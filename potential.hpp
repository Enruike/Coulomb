#pragma once
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<random>

//Clase Partícula
class Particle
{

    public:
    double pos_x, pos_y, pos_z; //Posiición en x, y, z.
    double valence; //Valencia. De modo que q = v * e.
    signed char specie; //Si es especie 1 o 2, incluso 3 si es un macroión.
    double radius; //Radio de la partícula
    //signed short int index; //Índice de la partícula
    bool infinite; //Homogénea o inhomogénea.

};

double vol_frac; //Volume fraction.
extern double delta_gr; //!!!(Check if we can vary this delta through parameters file)!!!
int dim_gr; //Dimension of the grid should be start from 0.                 
char *char_array;
double *radius_array;
double *valence_array;
double * very_initial_positions;
short int * cells;
short int * species_array;

extern double temp;
extern double eps_r;
extern int num_particles;
extern int species;
extern double r_1;
extern double r_2;
extern double val_1;
extern double val_2;
extern double *positions;
extern double box_len;
extern double diameter;
extern int max_time_steps;
extern int min_eq_steps;
extern int msd_steps;
extern int histo_steps;
extern int tau_steps;
extern double dt;
extern double diff_c_red;
extern double diff_c;
extern double diff_c_red2;
extern double macro_valence;
extern int macro_num;
extern bool file_pos_gen;

extern bool read_parameters();
extern bool read_file_atom_pos();
extern bool pos_macro_gen();
extern void periodic_distance(double xi, double yi, double zi,
    double& x_pos, double& y_pos, double& z_pos, double* positions);
extern double new_pos_function(int indx, double dt, double * indx_positions, double new_pos[3], short int cells[3]);
extern void histogram_hr_tau(int num_particles, double * positions, short int * species_array,
        double *** HR);
extern void calculate_rhor_gr(double *** RHOR, double *** GR, double *** HR, double tau, unsigned int * atoms_per_specie, double * bin_vol);
extern void write_gr_rhor(FILE * file, char * file_name, int tau, int species, double * XR, double *** gr_rhor);