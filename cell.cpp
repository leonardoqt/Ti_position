#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "class.h"

using namespace std;

void cell :: clean()
{}

void cell :: setup_cell(vec aa, vec bb, vec cc, int n_atom, ifstream& input)
{
	string flag_pos="RIMCOORD";
	string tmp;
	double r0 = 3.3;
	// element for SrTiO3
	num_element = 3;
	e_l = new string[3];
	e_l[0] = "Sr";
	e_l[1] = "Ti";
	e_l[2] = "O";
	// cell parameters
	latt[0] = aa;
	latt[1] = bb;
	latt[2] = cc;
	// all atoms
	num_atom = n_atom;
	a_l = new atom[num_atom];
	// read in all atom, and find how many Ti
	getline(input,tmp);
	while(tmp.find(flag_pos) == string::npos)
		getline(input,tmp);
	num_Ti = 0;
	for (int t1=0; t1<num_atom; t1++)
	{
		input>>a_l[t1].symbol>>a_l[t1].pos.x[0]>>a_l[t1].pos.x[1]>>a_l[t1].pos.x[2];	//cartisian coordinates
		getline(input, tmp);
		for (int t2=0; t2<num_element; t2++)
			if (a_l[t1].symbol == e_l[t2])
			{
				a_l[t1].sym = t2;
				break;
			}
		if (a_l[t1].sym == 1)	// 1 is Ti
			num_Ti++;
	}
	// get oxygen neighbor
	Ti_index = new int[num_Ti];
	oct_center = new vec[num_Ti];
	oct = new atom*[num_Ti];
	oct_basis = new vec*[num_Ti];
	for(int t1=0; t1<num_Ti; t1++)
	{
		oct[t1] = new atom[6];	// each Ti has 6 neighbor
		oct_basis[t1] = new vec[3];
	}
	for(int t1=0, t2=0; t1<num_atom; t1++)
		if (a_l[t1].sym == 1)
		{
			Ti_index[t2] = t1;
			t2++;
		}
	for(int t1=0; t1<num_Ti; t1++)
		for(int t_at = 0, t_o = 0; t_at<num_atom; t_at++)
		{
			if(a_l[t_at].sym == 2)
			{
				for(int dx=-1; dx<2; dx++)
				for(int dy=-1; dy<2; dy++)
				for(int dz=-1; dz<2; dz++)
				{
					if( (a_l[t_at].pos + latt[0]*dx + latt[1]*dy + latt[2]*dz - a_l[Ti_index[t1]].pos).norm()<r0 )
					{
						oct[t1][t_o] = a_l[t_at];
						oct[t1][t_o].pos = a_l[t_at].pos + latt[0]*dx + latt[1]*dy + latt[2]*dz;
						t_o++;
					}
				}
			}
		}
}

void cell :: get_local_coord(vec xx[3])
{
	double lim=0.85;	//~30deg
	vec tmp;
	// find oct center
	for(int t1=0; t1<num_Ti; t1++)
	{
		oct_center[t1].clean();
		for(int t2=0; t2<6; t2++)
			oct_center[t1] = oct_center[t1] + oct[t1][t2].pos;
		oct_center[t1] = oct_center[t1]*(1.0/6);
	}
	// find three axis
	for (int t1=0; t1<3; t1++)
		for (int t2=0; t2<num_Ti; t2++)
			for (int n1=0; n1<6; n1++)
			for (int n2=0; n2<6; n2++)
			{
				if (n1 != n2)
				{
					tmp = oct[t2][n2].pos - oct[t2][n1].pos;
					if( (tmp*xx[t1])/tmp.norm()/xx[t1].norm() > lim )	//forward
					{
						oct_basis[t2][t1] = tmp*(1/tmp.norm());
					}
					else if ( -(tmp*xx[t1])/tmp.norm()/xx[t1].norm() > lim )	//backward
					{
						oct_basis[t2][t1] = tmp*(-1/tmp.norm());
					}
				}
			}
}

void cell :: move_Ti(vec dx)
{
	for (int t1=0; t1<num_Ti; t1++)
	{
		a_l[Ti_index[t1]].pos = oct_center[t1] + oct_basis[t1][0] * dx.x[0] + oct_basis[t1][1] * dx.x[1] + oct_basis[t1][2] * dx.x[2];
	}
}

void cell :: print_all()
{
	for (int t1=0; t1<num_atom; t1++)
	{
		if (a_l[t1].sym == 2)	//O
			cout<<a_l[t1].symbol<<setw(17)<<setprecision(9)<<fixed<<a_l[t1].pos.x[0]<<setw(14)<<setprecision(9)<<fixed<<a_l[t1].pos.x[1]<<setw(14)<<setprecision(9)<<fixed<<a_l[t1].pos.x[2]<<endl;
		else
			cout<<a_l[t1].symbol<<setw(16)<<setprecision(9)<<fixed<<a_l[t1].pos.x[0]<<setw(14)<<setprecision(9)<<fixed<<a_l[t1].pos.x[1]<<setw(14)<<setprecision(9)<<fixed<<a_l[t1].pos.x[2]<<" 0 0 0"<<endl;
	}
}
