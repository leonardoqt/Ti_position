#include <iostream>
#include <fstream>
#include <string>
#include "class.h"

using namespace std;

int main()
{
	ifstream input;
	string input_name;
	int num_atom;
	vec aa,bb,cc;
	vec local_coord[3];
	vec move;
	cell sto;

	cin>>input_name;
	cin>>num_atom;
	cin>>aa.x[0]>>aa.x[1]>>aa.x[2];
	cin>>bb.x[0]>>bb.x[1]>>bb.x[2];
	cin>>cc.x[0]>>cc.x[1]>>cc.x[2];
	cin>>local_coord[0].x[0]>>local_coord[0].x[1]>>local_coord[0].x[2];
	cin>>local_coord[1].x[0]>>local_coord[1].x[1]>>local_coord[1].x[2];
	cin>>local_coord[2].x[0]>>local_coord[2].x[1]>>local_coord[2].x[2];
	cin>>move.x[0]>>move.x[1]>>move.x[2];
	input.open(input_name);
	sto.setup_cell(aa,bb,cc,num_atom,input);
	sto.get_local_coord(local_coord);
	sto.move_Ti(move);
	sto.print_all();

	return 0;
}
