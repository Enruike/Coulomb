#pragma once
#include<stdio.h>
#include<iostream>

/*Parameters*/
int particles; //número de partículas.
int species; //Número de especies.
double r_1; //Radio para la especie 1.
double r_2; //Radio para la especie 2.
double val_1; //Valencia para la especie 1.
double val_2; //Valencia para la especie 2.
double rp_d; //Repulsive core distance.
double diff_c; //Diffusion coefficient.
bool disol; //Homogénea o inhomogénea (Disolución infinita).
double box_len; //Box length.
double dt; //Delta time dt.
int max_eq_time; //Max time for equilibration.
int max_t_steps; //Max time steps.
int msd_steps; //Saving steps for Mean Square Displacement.
int gr_steps; //Saving steps for g(r).
double eps_r; //Relative epsilon.
double temp; //Temperature in Kelvin(K).

/*Parameters for position file*/
int dim; //Dimensionality.


/*Functions*/
bool read_parameters(); //Reading parameters function;
bool read_file_atom_pos(); //Reading the position file.