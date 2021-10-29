#include "vec3f.h"
#include "BBox.h"
#include "triFace.h"

using namespace std;

int readObj(string objPath, int& vNum, int& fNum, triFace* faces){
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
	}
	delete[] point;
	delete[] face;
	return 0;
}



int main(){
	int vNum, fNum;
	triFace* faces = new triFace;
	int res = readObj("flag-2000-changed.obj", vNum, fNum, faces);
	cout << vNum << " " << fNum << endl;
	if(res){
		cout << "read obj failed" << endl;
	}
	vec3f point11, point12, point13, point21, point22, point23;
	int collisionNum = 0;
	for(int i = 0; i < fNum; ++i){
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
	}
	delete[] faces->points1;
	delete[] faces->points2;
	delete[] faces->points3;
	delete[] faces->gravityPoints;
	delete faces;
	return 0;
}