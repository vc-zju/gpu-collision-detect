#include "BBox.h"
#include "vec3f.h"
using namespace std;

int main(){
    vec3f v1(0.5, 0.5, 0.5);
    vec3f v2(1.5, 1.5, 1.5);
    vec3f v3(1, 1, 1);
    vec3f v4(2, 2, 2);
    BBox box1(v1, v3);
    BBox box2(v2, v4);
    box1.merge(box2);
    cout << box1 << box2;
    return 0;
}