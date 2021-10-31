#pragma once
#include <algorithm>
#include "Node.h"
using namespace std;

struct sorted{
    unsigned int sortedMortonCode;
    int sortedObjectID;
};

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
unsigned int morton3D(double x, double y, double z)
{
    x = min(max(x * 1024.0, 0.0), 1023.0);
    y = min(max(y * 1024.0, 0.0), 1023.0);
    z = min(max(z * 1024.0, 0.0), 1023.0);
    unsigned int xx = expandBits((unsigned int)x);
    unsigned int yy = expandBits((unsigned int)y);
    unsigned int zz = expandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}

int clz(unsigned int code){
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

int findSplit( sorted*       sortedCodeAndID,
               int           first,
               int           last)
{
    // Identical Morton codes => split the range in the middle.

    unsigned int firstCode = sortedCodeAndID[first].sortedMortonCode;
    unsigned int lastCode = sortedCodeAndID[last].sortedMortonCode;

    if (firstCode == lastCode)
        return (first + last) >> 1;

    // Calculate the number of highest bits that are the same
    // for all objects, using the count-leading-zeros intrinsic.

    //int commonPrefix = __clz(firstCode ^ lastCode);
    int commonPrefix = clz(firstCode ^ lastCode);

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
            unsigned int splitCode = sortedCodeAndID[newSplit].sortedMortonCode;
            int splitPrefix = clz(firstCode ^ splitCode);
            if (splitPrefix > commonPrefix)
                split = newSplit; // accept proposal
        }
    }
    while (step > 1);

    return split;
}

Node* generateHierarchy( triFace*      faces,
                         sorted*       sortedCodeAndID,
                         int           first,
                         int           last)
{
    // Single object => create a leaf node.

    if (first == last)
        return new LeafNode(sortedCodeAndID[first].sortedObjectID, faces->box[sortedCodeAndID[first].sortedObjectID]);

    // Determine where to split the range.

    int split = findSplit(sortedCodeAndID, first, last);

    // Process the resulting sub-ranges recursively.

    Node* childA = generateHierarchy(faces, sortedCodeAndID,
                                     first, split);
    Node* childB = generateHierarchy(faces, sortedCodeAndID,
                                     split + 1, last);
    return new InternalNode(childA, childB, first, last);
}

/*Node* generateHierarchy( unsigned int* sortedMortonCodes,
                         int*          sortedObjectIDs,
                         int           numObjects)
{
    LeafNode* leafNodes = new LeafNode[numObjects];
    InternalNode* internalNodes = new InternalNode[numObjects - 1];

    // Construct leaf nodes.
    // Note: This step can be avoided by storing
    // the tree in a slightly different way.

    for (int idx = 0; idx < numObjects; idx++) // in parallel
        leafNodes[idx].objectID = sortedObjectIDs[idx];

    // Construct internal nodes.

    for (int idx = 0; idx < numObjects - 1; idx++) // in parallel
    {
        // Find out which range of objects the node corresponds to.
        // (This is where the magic happens!)

        int2 range = determineRange(sortedMortonCodes, numObjects, idx);
        int first = range.x;
        int last = range.y;

        // Determine where to split the range.

        int split = findSplit(sortedMortonCodes, first, last);

        // Select childA.

        Node* childA;
        if (split == first)
            childA = &leafNodes[split];
        else
            childA = &internalNodes[split];

        // Select childB.

        Node* childB;
        if (split + 1 == last)
            childB = &leafNodes[split + 1];
        else
            childB = &internalNodes[split + 1];

        // Record parent-child relationships.

        internalNodes[idx].childA = childA;
        internalNodes[idx].childB = childB;
        childA->parent = &internalNodes[idx];
        childB->parent = &internalNodes[idx];
    }

    // Node 0 is the root.

    return &internalNodes[0];
}*/

