#pragma once
#include "BBox_gpu.h"


struct Node{
    __host__ __device__ Node():leftChild(nullptr), rightChild(nullptr){}
    __host__ __device__ Node(Node* leftChild, Node* rightChild): leftChild(leftChild), rightChild(rightChild){}
    __host__ __device__ virtual bool isLeaf() = 0;
    __host__ __device__ virtual bool isInternal() = 0;
    __host__ __device__ virtual BBox* getNodeBox() = 0;
    __host__ __device__ virtual int getIndex() = 0;
    Node* leftChild;
    Node* rightChild;
};

struct LeafNode: public Node{
    __host__ __device__ LeafNode(int sortedObjectID, const BBox& box):index(sortedObjectID), box(box){}
    __host__ __device__ bool isLeaf() {return true;}
    __host__ __device__ bool isInternal() {return false;}
    __host__ __device__ BBox* getNodeBox() {return &box;}
    __host__ __device__ int getIndex() {return index;}
    int index;
    BBox box;
};

struct InternalNode: public Node{
    __host__ __device__ InternalNode(Node* leftChild, Node* rightChild, int start, int end):Node(leftChild, rightChild), start(start), end(end){
        box = box_merge(*(leftChild->getNodeBox()), *(rightChild->getNodeBox()));
    }
    __host__ __device__ bool isLeaf() {return false;}
    __host__ __device__ bool isInternal() {return true;}
    __host__ __device__ BBox* getNodeBox() {return &box;}
    __host__ __device__ int getIndex() {return -1;}
    int start, end;
    BBox box;
};