#pragma once
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<random>
#include<ctime>

/* Extern variables */

extern double* positions;
extern char* char_array;
extern double *radius_array;
extern double *valence_array;

extern int num_particles;
extern double half_box;
extern double box_len;
extern double Bjerrum_red;
extern double Bjerrum_len;
extern double diff_c;
extern double diff_c_red;
extern double diff_c_red2;
extern double d_rc; //Sigma for the rc potential.
extern double delta_gr;
extern int dim_gr;
extern int diag_grid;
extern int species;

//Repulsive Core
double energy_rc_i_all(int indx, int np_total);
double erc(double r, double ri, double rj);

//Electrostatic energy
double energy_el_i_all(int indx, int np_total);
double eew(double r, double ri, double rj);

//Distance calculator function
void periodic_distance(double xi, double yi, double zi,
    double& x_pos, double& y_pos, double& z_pos, double* positions);

//Box Muller Function
double box_muller(double num_1, double num_2);
//Movie generator function
void HOOMD_xml_generator();
//New position function
double new_pos_function(int indx, double dt, double * indx_positions, double new_pos[3], short int cells[3]);
double frc(double rij, double ri, double rj);
double few0(double rij, double ri, double rj);
// Random generator and Box-Muller function.
double random_muller();
//check movement
double mi_after_move(double& x, double& y, double& z, short int& cell_x, short int& cell_y, short int& cell_z);
//Histogram Function
void histogram_hr_tau(int num_particles, double * positions, short int * species_array,
        double *** HR);
void calculate_rhor_gr(double *** RHOR, double *** GR, double *** HR, double tau, unsigned int * atoms_per_specie, double * bin_vol);
void write_gr_rhor(FILE * file, char * file_name, int tau, int species, double * XR, double *** gr_rhor);

void pos_gen(double pos[3]);

void macro_histo_f(double * histo, double macro_pos[3], const double diagonal[3], const double mag);
void write_hist_macro_f(FILE * file, char * file_name, int tau, double * XR, double * histogram);