#pragma once
#include<stdio.h>
#include<iostream>
#include<random>
#include<ctime>

/*Parameters*/
int num_particles; //número de partículas.
int species; //Número de especies.
double r_1; //Radio para la especie 1.
double r_2; //Radio para la especie 2.
double val_1; //Valencia para la especie 1.
double val_2; //Valencia para la especie 2.
double rp_d; //Repulsive core distance.
double d_rc; //Sigma for the rc potential.
double diff_c; //Diffusion coefficient.
double diff_c_red;
double diff_c_red2;
bool disol; //Homogénea o inhomogénea (Disolución infinita).
double box_len; //Box length.                              !!!(Check if it is in reduced units)!!!
double half_box; //Half box length.
double dt; //Delta time dt.
int min_eq_steps; //Minimum time steps for equilibration.
int max_time_steps; //Max time steps.
int msd_steps; //Saving steps for Mean Square Displacement.
int histo_steps; //Steps for calculating histogram.
int tau_steps; //Steps for calculating g(r) and rho(r).
double eps_r; //Relative epsilon.
double temp; //Temperature in Kelvin(K).
double diameter; //                                             !!!(Check this diameter)!!!
double delta_gr; //Delta Grid.
double diag_grid; //Size of diagonal grid.
bool macro_flag;
int macro_num;
double macro_valence;
double macro_radius;
double macro_rc;
bool file_pos_gen;

//Spring constants.
double K_spring = 1e-11; // Spring Constant [N/m]
double K_red = K_spring / 1e-10; // Reduced Spring Constant [N/Angstrom]

//Dynamic array memory
double* positions;

/*Functions*/
bool read_parameters(); //Reading parameters function;
bool read_file_atom_pos(); //Reading the position file.
bool pos_macro_gen();
bool pos_macro_read();
extern void pos_gen(double [3]);