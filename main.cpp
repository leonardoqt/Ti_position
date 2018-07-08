#include <iostream>
#include <fstream>
#include <string>
#include "def.h"

using namespace std;

int main()
{
	ifstream input;
	string input_name;
	int ti, o0,o1,o2,o3,o4,o5;
	int num_atom, num_Ti;
	double x,y,z;
	int t1;
	cell sto;

	cin>>input_name;
	cin>>x>>y>>z;
	cin>>num_atom>>num_Ti;
	input.open(input_name);
	sto.get_lattice(x,y,z);
	sto.get_num_atom_Ti(num_atom,num_Ti);
	// find neighobr of Ti
	for (t1=0; t1< num_Ti; t1++)
	{
		cin>>ti>>o0>>o1>>o2>>o3>>o4>>o5;
		sto.label_Ti_O(ti,o0,o1,o2,o3,o4,o5);
	}
	// periodic correction for O around Ti
	while(!cin.eof())
	{
		cin>>ti>>o0>>x>>y>>z;
		sto.O_correction(ti,o0,x,y,z);
	}
	// read position
	for (t1=0; t1<7; t1++)
	{
		sto.read_coord(input);
		sto.update_oct();
		sto.project_Ti();
//		cout<<t1<<'\t'<<sto.Ti_proj[0].x+sto.Ti_proj[1].x+sto.Ti_proj[2].x+sto.Ti_proj[3].x<<'\t'<<sto.Ti_proj[0].y+sto.Ti_proj[1].y+sto.Ti_proj[2].y+sto.Ti_proj[3].y<<'\t'<<sto.Ti_proj[0].z+sto.Ti_proj[1].z+sto.Ti_proj[2].z+sto.Ti_proj[3].z<<endl;
		cout<<t1+1<<'\t'<<sto.position_Ti(2).x<<'\t'<<sto.position_Ti(2).y<<'\t'<<sto.position_Ti(2).z<<endl;
	}

//	sto.print_neighbor();
//	sto.print_O_correction();
//	sto.print_coord();
//	sto.print_ti_center();
//	sto.print_ti_basis();
//	sto.print_proj();
	return 0;
}
