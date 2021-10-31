#pragma once
#include "BBox.h"


struct Node{
    Node():leftChild(nullptr), rightChild(nullptr){}
    Node(Node* leftChild, Node* rightChild): leftChild(leftChild), rightChild(rightChild){}
    virtual bool isLeaf() = 0;
    virtual bool isInternal() = 0;
    virtual BBox* getNodeBox() = 0;
    virtual int getIndex() = 0;
    Node* leftChild;
    Node* rightChild;
};

struct LeafNode: public Node{
    LeafNode(int sortedObjectID, const BBox& box):index(sortedObjectID), box(box){}
    bool isLeaf() {return true;}
    bool isInternal() {return false;}
    BBox* getNodeBox() {return &box;}
    int getIndex() {return index;}
    int index;
    BBox box;
};

struct InternalNode: public Node{
    InternalNode(Node* leftChild, Node* rightChild, int start, int end):Node(leftChild, rightChild), start(start), end(end){
        box = box_merge(*(leftChild->getNodeBox()), *(rightChild->getNodeBox()));
    }
    bool isLeaf() {return false;}
    bool isInternal() {return true;}
    BBox* getNodeBox() {return &box;}
    int getIndex() {return -1;}
    int start, end;
    BBox box;
};