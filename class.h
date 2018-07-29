#include <cstring>
#include <fstream>
#ifndef MY_CLASS
#define MY_CLASS

class vec;
class atom;
class cell;

class vec
{
public:
	double x[3];
friend atom;
friend cell;
public:
	void import(double *);
	void clean();			// reset value to zero
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec & operator=(const vec&);
	vec & operator=(double*);	// can replace import
	double operator*(const vec&);
	double norm();
	//debug
	void print();
};

class atom
{
public:
	std::string symbol;
	int sym;
	vec pos;
friend cell;
public:
	atom & operator=(const atom&);
};

class cell
{
public:
	// atom related
	int num_atom;
	int num_Ti;
	int *Ti_index;
	int num_element;
	atom * a_l;			// atom list
	std::string * e_l;		// element list
	atom ** oct;	// neighbor list
	// cell related
	vec latt[3];
	// local coordinates
	vec *oct_center;
	vec **oct_basis;
public:
	// end of work
	void clean();
	// setupf of cell
	void setup_cell(vec aa, vec bb, vec cc, int n_atom, std::ifstream& input);
	void get_local_coord(vec xx[3]);
	// move atom
	void move_Ti(vec dx);
	// print
	void print_all();
};

#endif
