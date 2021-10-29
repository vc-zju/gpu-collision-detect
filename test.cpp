#include "BBox.h"
#include "vec3f.h"
#include "Node.h"
#include "triFace.h"
using namespace std;

int main(){
    vec3f v1(0.5, 0.5, 0.5);
    vec3f v2(1.5, 1.5, 1.5);
    vec3f v3(1, 1, 1);
    vec3f v4(2, 2, 2);
    BBox box1(v1, v2);
    BBox box2(v3, v4);
    cout << box_contact(box1, box2) << endl;
    return 0;
}