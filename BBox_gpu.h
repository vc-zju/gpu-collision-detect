#pragma once

#include <iostream>
#include "vec3f_gpu.h"
#include "triFace_gpu.h" 

struct BBox 
{ 
    vec3f min;
    vec3f max;
    __host__ __device__ BBox():min(-1, -1, -1), max(-1, -1, -1){}
    __host__ __device__ BBox(const vec3f& min, const vec3f& max): min(min), max(max){}
    __host__ __device__ BBox& merge(const BBox& box); 
};

__host__ __device__ inline BBox& BBox::merge(const BBox& box){
    vmin(this->min, box.min);
    vmax(this->max, box.max);
    return *this;
}

// inline std::ostream& operator<<(std::ostream&os, const BBox& box){
//     return std::cout << box.min << box.max;
// }

__host__ __device__ inline bool box_contact(const BBox& box1, const BBox& box2){
    if(box1.max.x > box2.min.x && box1.max.y > box2.min.y && box1.max.z > box2.min.z){
        if(box2.max.x > box1.min.x && box2.max.y > box1.min.y && box2.max.z > box1.min.z){
            return true;
        }
    }
    return false;
}

__host__ __device__ inline BBox box_merge(BBox box1, const BBox& box2){
    vmin(box1.min, box2.min);
    vmax(box1.max, box2.max);
    return box1;
}

__host__ __device__ BBox getBox(const triFace& face, const int& faceIndex){
    vec3f min ,max;
    min.x = fmin(face.points1[faceIndex].x, face.points2[faceIndex].x, face.points3[faceIndex].x);
    min.y = fmin(face.points1[faceIndex].y, face.points2[faceIndex].y, face.points3[faceIndex].y);
    min.z = fmin(face.points1[faceIndex].z, face.points2[faceIndex].z, face.points3[faceIndex].z);
    max.x = fmax(face.points1[faceIndex].x, face.points2[faceIndex].x, face.points3[faceIndex].x);
    max.y = fmax(face.points1[faceIndex].y, face.points2[faceIndex].y, face.points3[faceIndex].y);
    max.z = fmax(face.points1[faceIndex].z, face.points2[faceIndex].z, face.points3[faceIndex].z);
    return BBox(min, max);
};
