#include <iostream>
#ifndef CELL
#define CELL
class coord;
class cell;

class coord
{
public:
	double x,y,z;
};

class cell
{
private:
	coord lattice;		// lattice vector
	int num_atom;		// number of atoms
	coord *atom_f;		// atom coordiantes in fractional coordiates
	coord *atom;		// atom coordiantes
	int num_Ti;			// number of Ti
	int *Ti;			// array of Ti
	int **O;			// O neighbor of Ti
	coord **O_cor;		// O neighbor position correction
	coord **oct_basis;	// basis expanded by oct.
	coord *oct_center;	// center of oct.
	coord *Ti_proj;		// Ti position project to oct. basis
	coord *Ti_sphere;	// Ti position project to oct. in sphereical coordiante
	char **element;		// element symbol
public:
	cell();
	~cell();
	// setup of cell
	void get_lattice(double a, double b, double c);
	void get_num_atom_Ti(int numatom, int numTi);
	void label_Ti_O(int ti, int o0,int o1,int o2,int o3,int o4,int o5);
	void O_correction(int ti, int o, double x, double y, double z);
	// read in data
	void read_coord(std::ifstream &input);
	void update_oct();
	// project Ti position to basis
	void project_Ti();
	// calculate new Ti_position based on defined devation
	void new_Ti_position(double r, double theta, double phi);
	// export data
	coord position_Ti(int label);	// average position of Ti in oct. basis
	// debug
	void print_neighbor();
	void print_O_correction();
	void print_coord();
	void print_ti_center();
	void print_ti_basis();
	void print_proj();
};
#endif
