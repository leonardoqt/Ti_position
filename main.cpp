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
	for (int t1=0; t1<8; t1++)
	{
		for(int t2=0; t2<6; t2++)
			cout<<sto.oct[t1][t2].symbol<<'\t'<<sto.oct[t1][t2].pos.x[0]<<'\t'<<sto.oct[t1][t2].pos.x[1]<<'\t'<<sto.oct[t1][t2].pos.x[2]<<endl;
		cout<<sto.a_l[sto.Ti_index[t1]].symbol<<'\t'<<sto.a_l[sto.Ti_index[t1]].pos.x[0]<<'\t'<<sto.a_l[sto.Ti_index[t1]].pos.x[1]<<'\t'<<sto.a_l[sto.Ti_index[t1]].pos.x[2]<<endl;
	}
}
