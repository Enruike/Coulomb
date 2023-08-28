#pragma once
#include<math.h>
#include<stdio.h>

//Clase Partícula
class Particle
{

    public:
    double pos_x, pos_y, pos_z; //Posiición en x, y, z.
    double valence; //Valencia. De modo que q = v * e.
    signed char specie; //Si es especie 1 o 2, incluso 3 si es un macroión.
    double radius; //Radio de la partícula
    signed short int index; //Índice de la partícula
    bool infinite; //Homogénea o inhomogénea.

};

extern double temp;
extern double eps_r;

extern bool read_parameters();