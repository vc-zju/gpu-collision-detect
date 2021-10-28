#pragma once
#include "BBox.h"


struct Node{
    Node* leftChild;
    Node* rightChild;
    BBox box;
    bool isLeafFlag;
    int faceIndex;
    bool isLeaf() {return isLeafFlag;}
};