#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iomanip>
#include "def.h"

using namespace std;
//-----------------------------
// setupf of cell
cell :: cell()
{
	atom_f = nullptr;
	atom = nullptr;
	Ti = nullptr;
	O = nullptr;
	O_cor = nullptr;
	oct_basis = nullptr;
	oct_center = nullptr;
	Ti_proj = nullptr;
	Ti_sphere = nullptr;
	element = nullptr;
}

cell :: ~cell()
{
	if (atom_f != nullptr)
		delete[] atom_f;
	if (atom != nullptr)
		delete[] atom;
	if (Ti != nullptr)
		delete[] Ti;
	if (O != nullptr)
	{
		for (int t1=0; t1<num_Ti; t1++)
			if (O[t1] != nullptr)
				delete[] O[t1];
		delete[] O;
	}
	if (O_cor != nullptr)
	{
		for (int t1=0; t1<num_Ti; t1++)
			if (O_cor[t1] != nullptr)
				delete[] O_cor[t1];
		delete[] O_cor;
	}
	if (oct_basis != nullptr)
	{
		for (int t1=0; t1<num_Ti; t1++)
			if (oct_basis[t1] != nullptr)
				delete[] oct_basis[t1];
		delete[] oct_basis;
	}
	if (oct_center != nullptr)
		delete[] oct_center;
	if (Ti_proj != nullptr)
		delete[] Ti_proj;
	if (Ti_sphere != nullptr)
		delete[] Ti_sphere;
	if (element != nullptr)
	{
		for (int t1=0; t1<num_Ti; t1++)
			if (element[t1] != nullptr)
				delete[] element[t1];
		delete[] element;
	}
}

void cell :: get_lattice(double a, double b, double c)
{
	lattice.x = a;
	lattice.y = b;
	lattice.z = c;
}

void cell :: get_num_atom_Ti(int numatom, int numTi)
{
	num_atom = numatom;
	num_Ti = numTi;
	if (atom_f == nullptr)
		atom_f = new coord[num_atom];
	if (atom == nullptr)
		atom = new coord[num_atom];
	if (Ti == nullptr)
		Ti = new int[num_Ti];
	if (O == nullptr)
	{
		O = new int*[num_Ti];
		for (int t1=0; t1<num_Ti; t1++)
		{
			O[t1] = new int[6];
			Ti[t1] = -1;
		}
	}
	if (O_cor == nullptr)
	{
		O_cor = new coord*[num_Ti];
		for (int t1=0; t1<num_Ti; t1++)
		{
			O_cor[t1] = new coord[6];
			for (int t2=0; t2<6; t2++)
				O_cor[t1][t2].x = O_cor[t1][t2].y = O_cor[t1][t2].z = 0;
		}
	}
	if (oct_basis == nullptr)
	{
		oct_basis = new coord*[num_Ti];
		for (int t1=0; t1<num_Ti; t1++)
			oct_basis[t1] = new coord[3];
	}
	if (oct_center == nullptr)
		oct_center = new coord[num_Ti];
	if (Ti_proj == nullptr)
		Ti_proj = new coord[num_Ti];
	if (Ti_sphere == nullptr)
		Ti_sphere = new coord[num_Ti];
	if (element == nullptr)
	{
		element = new char*[num_atom];
		for (int t1=0; t1<num_atom; t1++)
			element[t1] = new char[3];
	}
}

void cell :: label_Ti_O(int ti, int o0,int o1,int o2,int o3,int o4,int o5)
{	// label ti is initialized as -1, 0->x 1->y 2->z 3->-z 4->-y 5->-x
	int lb=0;
	while (Ti[lb] > -1)
		lb++;
	if (lb >= num_Ti)
	{
		cout<<"Too many Ti atoms assigned!"<<endl;
		abort();
	}
	Ti[lb] = ti-1;
	O[lb][0] = o0-1;
	O[lb][1] = o1-1;
	O[lb][2] = o2-1;
	O[lb][3] = o3-1;
	O[lb][4] = o4-1;
	O[lb][5] = o5-1;
}

void cell :: O_correction(int ti, int o, double x, double y, double z)
{
	int lb_ti,lb_o;
	for(lb_ti=0; lb_ti<num_Ti; lb_ti++)
		if (Ti[lb_ti] == ti - 1)
			break;
	if (lb_ti >= num_Ti)
	{
		cout<<"Ti("<<ti<<") not in the Ti list when making correction to O"<<endl;
		abort();
	}
	for(lb_o=0; lb_o<6; lb_o++)
		if (O[lb_ti][lb_o] == o - 1)
			break;
	if (lb_o >= 6)
	{
		cout<<"O("<<o<<") not in the O list of Ti("<<ti<<") when making correction to O"<<endl;
		abort();
	}
	O_cor[lb_ti][lb_o].x = x;
	O_cor[lb_ti][lb_o].y = y;
	O_cor[lb_ti][lb_o].z = z;
}
//-----------------------------------------------
// read in data
void cell :: read_coord(ifstream &input)
{
	string temp;
	string label="PRIMCOORD";
	getline(input, temp);
	while(temp.find(label) == string::npos)
		getline(input, temp);
	for(int t1=0; t1<num_atom; t1++)
	{
		input>>element[t1]>>atom[t1].x>>atom[t1].y>>atom[t1].z;
		atom_f[t1].x = atom[t1].x/lattice.x;
		atom_f[t1].y = atom[t1].y/lattice.y;
		atom_f[t1].z = atom[t1].z/lattice.z;
	}
}

void cell :: update_oct()
{
	double norm;
	for (int t1=0; t1<num_Ti; t1++)
	{
		oct_basis[t1][0].x = ((atom_f[O[t1][0]].x+O_cor[t1][0].x)-(atom_f[O[t1][5]].x+O_cor[t1][5].x))*lattice.x;
		oct_basis[t1][0].y = ((atom_f[O[t1][0]].y+O_cor[t1][0].y)-(atom_f[O[t1][5]].y+O_cor[t1][5].y))*lattice.y;
		oct_basis[t1][0].z = ((atom_f[O[t1][0]].z+O_cor[t1][0].z)-(atom_f[O[t1][5]].z+O_cor[t1][5].z))*lattice.z;
		oct_basis[t1][1].x = ((atom_f[O[t1][1]].x+O_cor[t1][1].x)-(atom_f[O[t1][4]].x+O_cor[t1][4].x))*lattice.x;
		oct_basis[t1][1].y = ((atom_f[O[t1][1]].y+O_cor[t1][1].y)-(atom_f[O[t1][4]].y+O_cor[t1][4].y))*lattice.y;
		oct_basis[t1][1].z = ((atom_f[O[t1][1]].z+O_cor[t1][1].z)-(atom_f[O[t1][4]].z+O_cor[t1][4].z))*lattice.z;
		oct_basis[t1][2].x = ((atom_f[O[t1][2]].x+O_cor[t1][2].x)-(atom_f[O[t1][3]].x+O_cor[t1][3].x))*lattice.x;
		oct_basis[t1][2].y = ((atom_f[O[t1][2]].y+O_cor[t1][2].y)-(atom_f[O[t1][3]].y+O_cor[t1][3].y))*lattice.y;
		oct_basis[t1][2].z = ((atom_f[O[t1][2]].z+O_cor[t1][2].z)-(atom_f[O[t1][3]].z+O_cor[t1][3].z))*lattice.z;
		norm = sqrt(oct_basis[t1][0].x * oct_basis[t1][0].x + oct_basis[t1][0].y * oct_basis[t1][0].y + oct_basis[t1][0].z * oct_basis[t1][0].z);
		oct_basis[t1][0].x = oct_basis[t1][0].x / norm;
		oct_basis[t1][0].y = oct_basis[t1][0].y / norm;
		oct_basis[t1][0].z = oct_basis[t1][0].z / norm;
		norm = sqrt(oct_basis[t1][1].x * oct_basis[t1][1].x + oct_basis[t1][1].y * oct_basis[t1][1].y + oct_basis[t1][1].z * oct_basis[t1][1].z);
		oct_basis[t1][1].x = oct_basis[t1][1].x / norm;
		oct_basis[t1][1].y = oct_basis[t1][1].y / norm;
		oct_basis[t1][1].z = oct_basis[t1][1].z / norm;
		norm = sqrt(oct_basis[t1][2].x * oct_basis[t1][2].x + oct_basis[t1][2].y * oct_basis[t1][2].y + oct_basis[t1][2].z * oct_basis[t1][2].z);
		oct_basis[t1][2].x = oct_basis[t1][2].x / norm;
		oct_basis[t1][2].y = oct_basis[t1][2].y / norm;
		oct_basis[t1][2].z = oct_basis[t1][2].z / norm;

		oct_center[t1].x = ((atom_f[O[t1][0]].x+O_cor[t1][0].x)+(atom_f[O[t1][1]].x+O_cor[t1][1].x)+(atom_f[O[t1][2]].x+O_cor[t1][2].x)+(atom_f[O[t1][3]].x+O_cor[t1][3].x)+(atom_f[O[t1][4]].x+O_cor[t1][4].x)+(atom_f[O[t1][5]].x+O_cor[t1][5].x))*lattice.x / 6;
		oct_center[t1].y = ((atom_f[O[t1][0]].y+O_cor[t1][0].y)+(atom_f[O[t1][1]].y+O_cor[t1][1].y)+(atom_f[O[t1][2]].y+O_cor[t1][2].y)+(atom_f[O[t1][3]].y+O_cor[t1][3].y)+(atom_f[O[t1][4]].y+O_cor[t1][4].y)+(atom_f[O[t1][5]].y+O_cor[t1][5].y))*lattice.y / 6;
		oct_center[t1].z = ((atom_f[O[t1][0]].z+O_cor[t1][0].z)+(atom_f[O[t1][1]].z+O_cor[t1][1].z)+(atom_f[O[t1][2]].z+O_cor[t1][2].z)+(atom_f[O[t1][3]].z+O_cor[t1][3].z)+(atom_f[O[t1][4]].z+O_cor[t1][4].z)+(atom_f[O[t1][5]].z+O_cor[t1][5].z))*lattice.z / 6;
	}
}

void cell :: project_Ti()
{
	for (int t1=0; t1< num_Ti; t1++)
	{
		Ti_proj[t1].x = (atom[Ti[t1]].x-oct_center[t1].x)*oct_basis[t1][0].x+(atom[Ti[t1]].y-oct_center[t1].y)*oct_basis[t1][0].y+(atom[Ti[t1]].z-oct_center[t1].z)*oct_basis[t1][0].z;
		Ti_proj[t1].y = (atom[Ti[t1]].x-oct_center[t1].x)*oct_basis[t1][1].x+(atom[Ti[t1]].y-oct_center[t1].y)*oct_basis[t1][1].y+(atom[Ti[t1]].z-oct_center[t1].z)*oct_basis[t1][1].z;
		Ti_proj[t1].z = (atom[Ti[t1]].x-oct_center[t1].x)*oct_basis[t1][2].x+(atom[Ti[t1]].y-oct_center[t1].y)*oct_basis[t1][2].y+(atom[Ti[t1]].z-oct_center[t1].z)*oct_basis[t1][2].z;
	}
	// find position in spherical coordinate
	for (int t1=0; t1< num_Ti; t1++)
	{
		Ti_sphere[t1].x = sqrt(Ti_proj[t1].x*Ti_proj[t1].x + Ti_proj[t1].y*Ti_proj[t1].y + Ti_proj[t1].z*Ti_proj[t1].z);
		Ti_sphere[t1].y = acos(Ti_proj[t1].z / Ti_sphere[t1].x) * 180 / 3.1415926536;
		Ti_sphere[t1].z = acos(Ti_proj[t1].x / Ti_sphere[t1].x / sin(acos(Ti_proj[t1].z / Ti_sphere[t1].x))) * 180 / 3.1415926536;
		if (Ti_proj[t1].y / sin(acos(Ti_proj[t1].z / Ti_sphere[t1].x)) < 0)
			Ti_sphere[t1].z = -Ti_sphere[t1].z + 360;
	}
}

void cell :: new_Ti_position(double r, double theta, double phi)
{
	for (int t1=0; t1<num_Ti; t1++)
	{
		atom[Ti[t1]].x = oct_center[t1].x + r * (sin(theta) * cos(phi) * oct_basis[t1][0].x + sin(theta) * sin(phi) * oct_basis[t1][1].x + cos(theta) * oct_basis[t1][2].x);
		atom[Ti[t1]].y = oct_center[t1].y + r * (sin(theta) * cos(phi) * oct_basis[t1][0].y + sin(theta) * sin(phi) * oct_basis[t1][1].y + cos(theta) * oct_basis[t1][2].y);
		atom[Ti[t1]].z = oct_center[t1].z + r * (sin(theta) * cos(phi) * oct_basis[t1][0].z + sin(theta) * sin(phi) * oct_basis[t1][1].z + cos(theta) * oct_basis[t1][2].z);
	}
}

coord cell :: position_Ti(int label)
{
	coord pos;
	pos.x = pos.y = pos.z = 0;
	if (label == 1)
	{
		for(int t1=0; t1< num_Ti; t1++)
		{
			pos.x += Ti_proj[t1].x;
			pos.y += Ti_proj[t1].y;
			pos.z += Ti_proj[t1].z;
		}
	}
	else if (label == 2)
	{
		for(int t1=0; t1< num_Ti; t1++)
		{
			pos.x += Ti_sphere[t1].x;
			pos.y += Ti_sphere[t1].y;
			pos.z += Ti_sphere[t1].z;
		}
	}
	pos.x /= num_Ti;
	pos.y /= num_Ti;
	pos.z /= num_Ti;
	return pos;
}

//--------------------------
// debug
void cell :: print_neighbor()
{
	for (int t1=0; t1<num_Ti; t1++)
	{
		cout<<Ti[t1]+1;
		for (int t2=0; t2<6; t2++)
			cout<<' '<<O[t1][t2]+1;
		cout<<endl;
	}
}

void cell :: print_O_correction()
{
	for (int t1=0; t1<num_Ti; t1++)
		for (int t2=0; t2<6; t2++)
			if(O_cor[t1][t2].x !=0 || O_cor[t1][t2].y !=0 || O_cor[t1][t2].z !=0)
			{
				cout<<Ti[t1]+1<<' '<<O[t1][t2]+1<<' '<<O_cor[t1][t2].x<<' '<<O_cor[t1][t2].y<<' '<<O_cor[t1][t2].z<<endl;
			}
}

void cell :: print_coord()
{
	for (int t1=0; t1<num_atom; t1++)
		if (element[t1][1] == 0)
			cout<<element[t1]<<setw(18)<<fixed<<setprecision(9)<<atom[t1].x<<setw(15)<<setprecision(9)<<atom[t1].y<<setw(15)<<setprecision(9)<<atom[t1].z<<endl;
		else
			cout<<element[t1]<<setw(17)<<fixed<<setprecision(9)<<atom[t1].x<<setw(15)<<setprecision(9)<<atom[t1].y<<setw(15)<<setprecision(9)<<atom[t1].z<<endl;
}

void cell :: print_ti_center()
{
	for (int t1=0; t1<num_Ti; t1++)
	{
		cout<<atom[Ti[t1]].x<<' '<<atom[Ti[t1]].y<<' '<<atom[Ti[t1]].z<<'\t'<<oct_center[t1].x<<' '<<oct_center[t1].y<<' '<<oct_center[t1].z<<endl;
	}
}

void cell :: print_ti_basis()
{
	for (int t1=0; t1<num_Ti; t1++)
	{
		cout<<oct_basis[t1][0].x<<' '<<oct_basis[t1][0].y<<' '<<oct_basis[t1][0].z<<'\t';
		cout<<oct_basis[t1][1].x<<' '<<oct_basis[t1][1].y<<' '<<oct_basis[t1][1].z<<'\t';
		cout<<oct_basis[t1][2].x<<' '<<oct_basis[t1][2].y<<' '<<oct_basis[t1][2].z<<endl;
	}
}

void cell :: print_proj()
{
	for (int t1=0; t1<num_Ti; t1++)
		cout<<Ti_proj[t1].x<<'\t'<<Ti_proj[t1].y<<'\t'<<Ti_proj[t1].z<<endl;
}
