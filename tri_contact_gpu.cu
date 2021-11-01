#include <stdlib.h>
#include <time.h>
#include <stdio.h>

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

Node* buildBVH(triFace* faces, const int& fNum, const vec3f& min, const vec3f& max){
	vec3f extend = max - min;
	sorted* sortedCodeAndID = new sorted[fNum];
	for(int i = 0; i < fNum; ++i){
		sortedCodeAndID[i].sortedMortonCode = morton3D((faces->gravityPoints[i].x - min.x) / extend.x,\
													   (faces->gravityPoints[i].y - min.y) / extend.y,\
													   (faces->gravityPoints[i].z - min.z) / extend.z);
		sortedCodeAndID[i].sortedObjectID = i;
	}
	qsort(sortedCodeAndID, fNum, sizeof(sorted), compare);
	Node* root = generateHierarchy(faces, sortedCodeAndID, 0, fNum - 1);
	delete[] sortedCodeAndID;
	return root;
}

__device__ void traverseBVHGPU(vec3f* points1, vec3f* points2, vec3f* points3, BBox* bbox, const int& queryNum, Node* root){
	const BBox& box = bbox[queryNum];
	if(box_contact(box, *(root->getNodeBox()))){
		if(root->isLeaf()){
			int index = root->getIndex();
			if(queryNum < index && tri_contact(points1[queryNum], points2[queryNum], points3[queryNum], 
						   points1[index], points2[index], points3[index])){
				printf("%d %d\n", queryNum, index);
			}
		}
		else{
			Node* leftChild = root->leftChild;
			Node* rightChild = root->rightChild;
			traverseBVHGPU(points1, points2, points3, bbox, queryNum, leftChild);
			traverseBVHGPU(points1, points2, points3, bbox, queryNum, rightChild);
		}
	}
}

__global__ void traverseBVH(triFace* faces, Node* root, int fNum){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i < fNum){
		traverseBVHGPU(faces->points1, faces->points2, faces->points3, faces->box, i, root);
	}
}

__global__ void allocateOnDevice(triFace** faces){
	
}

__global__ void freeOnDevice(triFace* faces, Node** root){
	
}

int main(){
	int vNum, fNum;
	clock_t start, end;
	vec3f min(100, 100, 100);
	vec3f max(-100, -100, -100);
	triFace* faces = new triFace;
	Node* root;
	start = clock();
	int res = readObj("flag-2000-changed.obj", vNum, fNum, faces, min, max);
	if(res){
		cout << "read obj failed" << endl;
	}
	end = clock();
	printf("read time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	triFace* facesGPU = new triFace;
	// CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU), sizeof(triFace)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points1), fNum  * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points2), fNum  * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->points3), fNum  * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->gravityPoints), fNum  * sizeof(vec3f)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMalloc((void**)&(facesGPU->box), fNum  * sizeof(BBox)), "cudaMalloc failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points1, faces->points1, fNum  * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points1 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points2, faces->points2, fNum  * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points2 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->points3, faces->points3, fNum  * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy points3 failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->gravityPoints, faces->gravityPoints, fNum  * sizeof(vec3f), cudaMemcpyHostToDevice), "cudaMemcpy gravityPoints failed!\n", -2);
	CUDA_CHECK_CALL(cudaMemcpy(facesGPU->box, faces->box, fNum  * sizeof(BBox), cudaMemcpyHostToDevice), "cudaMemcpy box failed!\n", -2);
	start = clock();
	root = buildBVH(faces, fNum, min, max);
	end = clock();
	printf("build time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	start = clock();
	traverseBVH<<<(fNum + 1023) / 1024, 1024>>>(facesGPU, root, fNum);
	end = clock();
	printf("traverse time=%fms\n",(double)(end-start) * 1000 /CLK_TCK);
	cudaFree(facesGPU->points1);
	cudaFree(facesGPU->points2);
	cudaFree(facesGPU->points3);
	cudaFree(facesGPU->gravityPoints);
	cudaFree(facesGPU->box);
	delete facesGPU;
	/*vec3f point11, point12, point13, point21, point22, point23;
	int collisionNum = 0;
	for(int i = 120914; i < fNum; ++i){
		for(int j = i + 1; j < fNum; ++j){
			point11 = faces->points1[i];
			point12 = faces->points2[i];
			point13 = faces->points3[i];
			point21 = faces->points1[j];
			point22 = faces->points2[j];
			point23 = faces->points3[j];
			if(tri_contact(point11, point12, point13, point21, point22, point23)){
				if(!(point11.equal_abs(point21) || point11.equal_abs(point22) || point11.equal_abs(point23) ||\
				     point12.equal_abs(point21) || point12.equal_abs(point22) || point12.equal_abs(point23) ||\
					 point13.equal_abs(point21) || point13.equal_abs(point22) || point13.equal_abs(point23))){
					++collisionNum;
					cout << "#self contact found at (" << i << "," << j << ")" << endl;
					cout << point11 << endl;
					cout << point12 << endl;
					cout << point13 << endl;
					cout << point21 << endl;
					cout << point22 << endl;
					cout << point23 << endl;
				}	
			}
		}
	}*/
	delete[] faces->points1;
	delete[] faces->points2;
	delete[] faces->points3;
	delete[] faces->gravityPoints;
	delete[] faces->box;
	delete faces;
	return 0;
}