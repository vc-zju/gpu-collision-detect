#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>

#include "vec3f_gpu.h"
#include "BBox_gpu.h"
#include "triFace_gpu.h"
#include "Node_gpu.h"
#include "build_gpu.h"

using namespace std;

cudaError_t gCudaStatus;

#define CUDA_CHECK_CALL(fun, err_msg, return_code)					\
	gCudaStatus = fun;												\
	if(gCudaStatus != cudaSuccess){									\
		printf("error_code%d: %s", gCudaStatus, err_msg);			\
		return return_code;											\
	}

int readObj(string objPath, int& vNum, int& fNum, triFace* faces, vec3f& min, vec3f& max){
	ifstream infile(objPath);
	if(!infile.is_open()){
		return -1;
	}
	string line;
	vNum = 0;
	fNum = 0;
	while(getline(infile, line)){
		if(line[0] == 'v' && line[1] == ' '){
			++vNum; 
		}
		else if(line[0] == 'f' && line[1] == ' '){
			++fNum;
		}
	}
	vec3f* point = new vec3f[vNum];
	vec3f* face = new vec3f[fNum];
	faces->points1 = new vec3f[fNum];
	faces->points2 = new vec3f[fNum];
	faces->points3 = new vec3f[fNum];
	faces->gravityPoints = new vec3f[fNum];
	faces->box = new BBox[fNum];
	infile.close();
	infile.open(objPath);
	if(!infile.is_open()){
		return -1;
	}
	int vNumIndex = 0, fNUmIndex = 0;
	string str;
	double d[3];
	while(getline(infile, line)){
		if(line[0] == 'v' && line[1] == ' '){
			istringstream instr(line);
			instr >> str >> d[0] >> d[1] >> d[2];
			// cout << str << " " << d1 << " " << d2 << " " << d3 << endl;
			point[vNumIndex].set_value(d[0], d[1], d[2]);
			++vNumIndex; 
			if(d[0] > max.x){
				max.x = d[0];
			}
			if(d[0] < min.x){
				min.x = d[0];
			}
			if(d[1] > max.y){
				max.y = d[1];
			}
			if(d[1] < min.y){
				min.y = d[1];
			}
			if(d[2] > max.z){
				max.z = d[2];
			}
			if(d[2] < min.z){
				min.z = d[2];
			}
		}
		else if(line[0] == 'f' && line[1] == ' '){
			istringstream instr(line);
			instr >> str;
			for(int i = 0; i < 3; ++i){
				instr >> str;
				d[i] = atof(str.c_str());
			}
			face[fNUmIndex].set_value(d[0], d[1], d[2]);
			++fNUmIndex;
		}
	}
	infile.close();
	for(int i = 0; i < fNum; ++i){
		faces->points1[i] = point[static_cast<int>(face[i].x) - 1];
		faces->points2[i] = point[static_cast<int>(face[i].y) - 1];
		faces->points3[i] = point[static_cast<int>(face[i].z) - 1];
		faces->gravityPoints[i] = (faces->points1[i] + faces->points2[i] + faces->points3[i]) / 3;
		faces->box[i] = getBox(*faces, i);
	}
	delete[] point;
	delete[] face;
	return 0;
}

int compare(const void* a, const void* b){
	return ((sorted*)a)->sortedMortonCode - ((sorted*)b)->sortedMortonCode;
}

__global__ void calclulateMorton3D(vec3f* gravityPoints, int fNum, unsigned int* sortedMortonCode, int* sortedObjectID, const vec3f min, const vec3f extend){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < fNum){
		sortedMortonCode[i] = morton3D((gravityPoints[i].x - min.x) / extend.x, (gravityPoints[i].y - min.y) / extend.y, (gravityPoints[i].z - min.z) / extend.z);
		sortedObjectID[i] = i;
	}
}
 
void buildBVH(triFace* faces, triFace* facesGPU, const int& fNum, const vec3f& min, const vec3f& max, Node** root){
	vec3f extend = max - min;
	// sorted* sortedCodeAndID = new sorted[fNum];
	// sorted* sortedCodeAndIDGPU;
	// cudaMalloc((void**)&(sortedCodeAndIDGPU), fNum * sizeof(sorted));
	unsigned int* sortedMortonCode;
    int* sortedObjectID;
	cudaMalloc((void**)&(sortedMortonCode), fNum * sizeof(unsigned int));
	cudaMalloc((void**)&(sortedObjectID), fNum * sizeof(int));
	calclulateMorton3D<<<(fNum + 1023) / 1024, 1024>>>(facesGPU->gravityPoints, fNum, sortedMortonCode, sortedObjectID, min, extend);
	cudaThreadSynchronize();
	// cudaMemcpy(sortedCodeAndID, sortedCodeAndIDGPU, fNum * sizeof(sorted), cudaMemcpyDeviceToHost);
	thrust::device_ptr<unsigned int> dev_key_ptr(sortedMortonCode);
	thrust::device_ptr<int> dev_data_ptr(sortedObjectID);
	thrust::sort_by_key(dev_key_ptr, dev_key_ptr + fNum, dev_data_ptr);
	sortedMortonCode = thrust::raw_pointer_cast(dev_key_ptr);
	sortedObjectID = thrust::raw_pointer_cast(dev_data_ptr);
	// qsort(sortedCodeAndID, fNum, sizeof(sorted), compare);
	// cudaMemcpy(sortedCodeAndIDGPU, sortedCodeAndID, fNum * sizeof(sorted), cudaMemcpyHostToDevice);
	generateHierarchy(facesGPU->box, sortedMortonCode, sortedObjectID, fNum, root);
	// cudaFree(sortedCodeAndIDGPU);
	cudaFree(sortedMortonCode);
	cudaFree(sortedObjectID);
}

__device__ void traverseBVHGPU(vec3f* points1, vec3f* points2, vec3f* points3, BBox* bbox, const int& queryNum, Node* rootPtr){
	const BBox& box = bbox[queryNum];
	if(box_contact(box, *(rootPtr->getNodeBox()))){
		if(rootPtr->isLeaf()){
			int index = rootPtr->getIndex();
			if(queryNum < index && tri_contact(points1[queryNum], points2[queryNum], points3[queryNum], 
						   points1[index], points2[index], points3[index])){
				// printf("%d %d\n", queryNum, index);
			}
		}
		else{
			Node* leftChild = rootPtr->leftChild;
			Node* rightChild = rootPtr->rightChild;
			traverseBVHGPU(points1, points2, points3, bbox, queryNum, leftChild);
			traverseBVHGPU(points1, points2, points3, bbox, queryNum, rightChild);
		}
	}
}

__device__ void traverseIterativeBVHGPU( vec3f* points1, vec3f* points2, vec3f* points3, BBox* bbox, const int& queryNum, Node* root)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.
    Node* stack[80];
    Node** stackPtr = stack;
    *stackPtr++ = NULL; // push

    // Traverse nodes starting from the root.
    do
    {
        // Check each child node for overlap.
        Node* childL = root->leftChild;
        Node* childR = root->rightChild;
        bool overlapL = ( box_contact(bbox[queryNum], *(childL->getNodeBox())));
        bool overlapR = ( box_contact(bbox[queryNum], *(childR->getNodeBox())));

        // Query overlaps a leaf node => report collision.
        if (overlapL && childL->isLeaf()){
			int index = childL->getIndex();
			if(queryNum < index && tri_contact(points1[queryNum], points2[queryNum], points3[queryNum], 
						   points1[index], points2[index], points3[index])){
				printf("%d %d\n", queryNum, childL->getIndex());
			}
		}

        if (overlapR && childR->isLeaf()){
			int index = childR->getIndex();

			if(queryNum < index && tri_contact(points1[queryNum], points2[queryNum], points3[queryNum], 
						   points1[index], points2[index], points3[index])){
				printf("%d %d\n", queryNum, childR->getIndex());
			}
		}
        // Query overlaps an internal node => traverse.
        bool traverseL = (overlapL && !childL->isLeaf());
        bool traverseR = (overlapR && !childR->isLeaf());

        if (!traverseL && !traverseR)
            root = *--stackPtr; // pop
        else
        {
            root = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    }
    while (root != NULL);
}

__global__ void traverseBVH(vec3f* points1, vec3f* points2, vec3f* points3, BBox* box, Node** root, int fNum){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	Node* rootPrt = *root;
	if(i < fNum){
		// traverseBVHGPU(points1, points2, points3, box, i, rootPrt);
		traverseIterativeBVHGPU(points1, points2, points3, box, i, rootPrt);
	}
}


__global__ void freeOnDevice(Node** root){
	delete (*root);
	*root = nullptr;
}

int main(){
	size_t limit = 65536;//4096;	//set cuda stack size to avoid stackoverflow
	cudaDeviceSetLimit(cudaLimitStackSize, limit);
	cudaDeviceGetLimit(&limit, cudaLimitStackSize);
	printf("cudaLimitStackSize: %llu\n", limit);
	limit = 536870912;
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, limit);
	cudaDeviceGetLimit(&limit, cudaLimitMallocHeapSize);
	printf("cudaLimitMallocHeapSize: %llu\n", limit);
	int vNum, fNum;
	clock_t start, end;
	vec3f min(100, 100, 100);
	vec3f max(-100, -100, -100);
	triFace* faces = new triFace;
	Node** root;
	cudaMalloc((void **)&root, sizeof(Node*));
	start = clock();
	int res = readObj("flag-2000-changed.obj", vNum, fNum, faces, min, max);
	if(res){
		cout << "read obj failed" << endl;
	}
	end = clock();
	printf("read time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	triFace* facesGPU = new triFace;
	// CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU), sizeof(triFace)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points1), fNum * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points2), fNum * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points3), fNum * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->gravityPoints), fNum * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->box), fNum * sizeof(BBox)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points1, faces->points1, fNum * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points1 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points2, faces->points2, fNum * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points2 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points3, faces->points3, fNum * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points3 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->gravityPoints, faces->gravityPoints, fNum * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy gravityPoints failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->box, faces->box, fNum * sizeof(BBox), cudaMemcpyHostToDevice), "cudaMemcpy box failed!\n", -2);
	start = clock();
	buildBVH(faces, facesGPU, fNum, min, max, root);
	end = clock();
	printf("build time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	start = clock();
	traverseBVH<<<(fNum + 511) / 512, 512>>>(facesGPU->points1, facesGPU->points2, facesGPU->points3, facesGPU->box, root, fNum);
	cudaDeviceSynchronize();
	end = clock();
	printf("traverse time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	cudaFree(facesGPU->points1);
	cudaFree(facesGPU->points2);
	cudaFree(facesGPU->points3);
	cudaFree(facesGPU->gravityPoints);
	cudaFree(facesGPU->box);
	delete facesGPU;
	freeOnDevice<<<1, 1>>>(root);
	cudaFree(root);
	delete[] faces->points1;
	delete[] faces->points2;
	delete[] faces->points3;
	delete[] faces->gravityPoints;
	delete[] faces->box;
	delete faces;
	return 0;
}