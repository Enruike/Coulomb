#pragma once
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<random>

//Clase Partícula
class Particule
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
double delta_gr; //                     !!!(Check if we can vary this delta through parameters file)!!!
char *char_array;
double *radius_array;
double *valence_array;

extern double temp;
extern double eps_r;
extern int num_particules;
extern int species;
extern double r_1;
extern double r_2;
extern double val_1;
extern double val_2;
extern double *positions;
extern double box_len;
extern double diameter;

extern bool read_parameters();
extern bool read_file_atom_pos();