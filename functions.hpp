#pragma once
#include<stdlib.h>
#include<math.h>
#include<iostream>

/* Extern variables */

extern double* positions;
extern char* char_array;
extern double *radius_array;
extern double *valence_array;

extern int num_particules;
extern double half_box;
extern double box_len;
extern double Bjerrum_red;

//Repulsive Core
double energy_rc_i_all(int indx, int np_total);
double erc(double r, double ri, double rj);

//Electrostatic energy
double energy_el_i_all(int indx, int np_total);
double eew(double r, double ri, double rj);

//Box Muller Function
double box_muller(double num_1, double num_2);
//Movie generator function
void HOOMD_xml_generator();
