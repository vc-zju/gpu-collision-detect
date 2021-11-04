#pragma once
#include <algorithm>
#include "Node_gpu.h"
using namespace std;

struct sorted{
    unsigned int sortedMortonCode;
    int sortedObjectID;
};

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
__host__ __device__ unsigned int expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
__host__ __device__ unsigned int morton3D(double x, double y, double z)
{
    x = min(max(x * 1024.0, 0.0), 1023.0);
    y = min(max(y * 1024.0, 0.0), 1023.0);
    z = min(max(z * 1024.0, 0.0), 1023.0);
    unsigned int xx = expandBits((unsigned int)x);
    unsigned int yy = expandBits((unsigned int)y);
    unsigned int zz = expandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}

__host__ __device__ int clz(unsigned int code){
    int zeroNum = 0;
    unsigned int pattern = 0x80000000;
    for(int i = 0; i < 32; ++i){
        if((code & pattern) == 0){
            ++zeroNum;
            pattern = pattern >> 1;
        }
        else{
            return zeroNum;
        }
    }
    return zeroNum;
}

__device__ int findSplit( unsigned int* sortedMortonCode,
               int           first,
               int           last)
{
    // Identical Morton codes => split the range in the middle.

    unsigned int firstCode = sortedMortonCode[first];
    unsigned int lastCode = sortedMortonCode[last];

    if (firstCode == lastCode)
        return (first + last) >> 1;

    // Calculate the number of highest bits that are the same
    // for all objects, using the count-leading-zeros intrinsic.

    //int commonPrefix = __clz(firstCode ^ lastCode);
    int commonPrefix = __clz(firstCode ^ lastCode);

    // Use binary search to find where the next bit differs.
    // Specifically, we are looking for the highest object that
    // shares more than commonPrefix bits with the first one.

    int split = first; // initial guess
    int step = last - first;

    do
    {
        step = (step + 1) >> 1; // exponential decrease
        int newSplit = split + step; // proposed new position

        if (newSplit < last)
        {
            unsigned int splitCode = sortedMortonCode[newSplit];
            int splitPrefix = __clz(firstCode ^ splitCode);
            if (splitPrefix > commonPrefix)
                split = newSplit; // accept proposal
        }
    }
    while (step > 1);

    return split;
}

/*__host__ __device__ Node* generateHierarchy( BBox*      box,
                         sorted*       sortedCodeAndID,
                         int           first,
                         int           last)
{
    // Single object => create a leaf node.
    if (first == last){
        return new LeafNode(sortedCodeAndID[first].sortedObjectID, box[sortedCodeAndID[first].sortedObjectID]);
    }
        

    // Determine where to split the range.

    int split = findSplit(sortedCodeAndID, first, last);

    // Process the resulting sub-ranges recursively.

    Node* childA = generateHierarchy(box, sortedCodeAndID,
                                     first, split);
    Node* childB = generateHierarchy(box, sortedCodeAndID,
                                     split + 1, last);
    return new InternalNode(childA, childB, first, last);
}*/

__host__ __device__ inline int sign(int num){
    if(num >= 0){
        return 1;
    }
    else if(num < 0){
        return -1;
    }
}

__device__ inline int clzMorton(unsigned int* sortedMortonCode, int numObjects, int idx1, int idx2){
    if((idx2 >= 0) && (idx2 <= numObjects - 1)){
        if(((sortedMortonCode[idx1]) ^ (sortedMortonCode[idx2])) == 0){
            return (32 + __clz(idx1 ^ idx2));
        }
        else
            return __clz((sortedMortonCode[idx1]) ^ (sortedMortonCode[idx2]));
    }  
    else
        return -1;
}

__device__ void determineRange(unsigned int* sortedMortonCode, int numObjects, int idx, int* first, int* last){
    int d = sign(clzMorton(sortedMortonCode, numObjects, idx, idx + 1) - clzMorton(sortedMortonCode, numObjects, idx, idx - 1));
    int dMin = clzMorton(sortedMortonCode, numObjects, idx, idx - d);
    int bound = 2;
    while(clzMorton(sortedMortonCode, numObjects, idx, idx + d * bound) > dMin){
        bound *= 2;
    }
    int boundAnother = 0;
    for(int t = bound / 2; t >= 1; t /= 2){
        if(clzMorton(sortedMortonCode, numObjects, idx, idx + (boundAnother + t) * d) > dMin){
            boundAnother = boundAnother + t;
        }
    }
    int j = idx + boundAnother * d;
    if(d > 0){
        *last = j;
        *first = idx;
    }
    else{
        *last = idx;
        *first = j;
    }
}

__global__ void allocateOnDevice(LeafNode** leafNodes, InternalNode** internalNodes, int numObjects, Node** root){
	*leafNodes = new LeafNode[numObjects];
    *internalNodes = new InternalNode[numObjects - 1];
    *root = &(*internalNodes)[0];
}

__global__ void assignLeafNodes(LeafNode** leafNodes, int* sortedObjectID, BBox* box,int numObjects){  
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    LeafNode* leafNodesPtr = *leafNodes;
    if(idx < numObjects){
        leafNodesPtr[idx].index = sortedObjectID[idx];
        leafNodesPtr[idx].box = box[sortedObjectID[idx]];
    }
}

__global__ void assignInternalNodes(LeafNode** leafNodes, InternalNode** internalNodes, unsigned int* sortedMortonCode, int* sortedObjectID, int numObjects){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    LeafNode* leafNodesPtr = *leafNodes;
    InternalNode* internalNodesPtr = *internalNodes;
    if(idx < numObjects - 1){
        // Find out which range of objects the node corresponds to.
        // (This is where the magic happens!)
        int first, last;
        determineRange(sortedMortonCode, numObjects, idx, &first, &last);
        // Determine where to split the range.

        int split = findSplit(sortedMortonCode, first, last);

        // Select childA.

        Node* childA;
        if (split == first)
            childA = &(leafNodesPtr[split]);
        else
            childA = &(internalNodesPtr[split]);

        // Select childB.

        Node* childB;
        if (split + 1 == last){
            childB = &(leafNodesPtr[split + 1]);
        }
        else
            childB = &(internalNodesPtr[split + 1]);

        // Record parent-child relationships.

        internalNodesPtr[idx].leftChild = childA;
        internalNodesPtr[idx].rightChild = childB;
        internalNodesPtr[idx].box = leafNodesPtr[first].box; 
        for(int i = first + 1; i <= last; ++i){
            internalNodesPtr[idx].box = box_merge(internalNodesPtr[idx].box, leafNodesPtr[i].box);
        }
    }
}

__global__ void debug(Node** root){
    printf("%.6f", (*root)->getNodeBox()->max.x);
    printf("%.6f", (*root)->getNodeBox()->max.y);
    printf("%.6f", (*root)->getNodeBox()->max.z);
    printf("%.6f", (*root)->getNodeBox()->min.x);
    printf("%.6f", (*root)->getNodeBox()->min.y);
    printf("%.6f", (*root)->getNodeBox()->min.z);
}

__host__ void generateHierarchy( BBox*      box,
                         unsigned int* sortedMortonCode, 
                         int*          sortedObjectID,
                         int           numObjects,
                         Node**        root)
{
    LeafNode** leafNodes;
    InternalNode** internalNodes;
    cudaMalloc((void **)&leafNodes, sizeof(LeafNode*));
    cudaMalloc((void **)&internalNodes, sizeof(InternalNode*));
    allocateOnDevice<<<1, 1>>>(leafNodes, internalNodes, numObjects, root);
    // Construct leaf nodes.
    // Note: This step can be avoided by storing
    // the tree in a slightly different way.
    assignLeafNodes<<<(numObjects + 1023) / 1024, 1024>>>(leafNodes, sortedObjectID, box, numObjects);   
    // Construct internal nodes.
    assignInternalNodes<<<(numObjects + 1023) / 1024, 1024>>>(leafNodes, internalNodes, sortedMortonCode, sortedObjectID, numObjects);
    cudaDeviceSynchronize();
}


