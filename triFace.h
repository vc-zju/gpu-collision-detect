#pragma once
#include "vec3f.h"

struct BBox;

struct triFace 
{
	vec3f* points1;
	vec3f* points2;
	vec3f* points3;
    vec3f* gravityPoints;
	BBox* box;
};