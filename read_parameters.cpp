#include "read_parameters.hpp"

bool read_parameters(){

    FILE* param = fopen("param.in", "r");

	if (param == (FILE*)NULL) {
		return false;
	}
    else{
        fscanf(param, "**** PARAMETERS FILE ****\n\n");
        fscanf(param, "Number of particles = %d\n", &particles);
        fscanf(param, "Number of species = %d\n", &species);
        fscanf(param, "Radius of species 1 = %lf\n", &r_1);
        fscanf(param, "Radius of species 2 = %lf\n", &r_2);
        fscanf(param, "Valence of species 1 = %lf\n", &val_1);
        fscanf(param, "Valence of species 2 = %lf\n", &val_2);
        fscanf(param, "Repulsive core distance = %lf\n", &rp_d);
        fscanf(param, "Diffusion Coefficient = %lf\n", &diff_c);
        fscanf(param, "Length of the box = %lf #Length in Angstroms\n", &box_len);
        fscanf(param, "Delta t = %lf #Time step in seconds\n", &dt);
        fscanf(param, "Max equilibration time = %d #Max time for equilibration\n", &max_eq_time);
        fscanf(param, "Max time steps = %d #Max number of steps\n", &max_t_steps);
        fscanf(param, "MSD steps = %d #Save MSD file at certain steps\n", &msd_steps);
        fscanf(param, "g(r) steps = %d #Saving steps for g(r)\n", &gr_steps);
        fscanf(param, "Infinite disolution = %d #0:Homogeneous or 1:Inhomogeneous\n", &disol);
        fscanf(param, "Epsilon r = %lf #Epsilon for electrolyte\n", &eps_r);
        fscanf(param, "Temperature = %lf #Temperature in Kelvin(K)\n", &temp);

    }
    
    char* flag;

    if(disol == 0){
        flag = "False";
    }
    else{
        flag = "True";
    }

    printf("**** PARAMETERS FILE ****\n\n\
    Number of particles = %d\n\
    Number of species = %d\n\
    Radius of species 1 = %.2lf\u00c5\n\
    Radius of species 2 = %.2lf\u00c5\n\
    Valence of species 1 = %.2lf\n\
    Valence of species 2 = %.2lf\n\
    Repulsive core distance = %.2lf\u00c5\n\
    Diffusion Coefficient = %.2e\n\
    Length of the box = %.2lf\u00c5\n\
    Delta t = %.2e s\n\
    Max equilibration time = %d\n\
    Max time steps = %d\n\
    MSD saving steps = %d\n\
    g(r) saving steps = %d\n\
    Infinite disolution = %s\n\
    Epsilon r = %.2lf\n\
    Temperature = %.2lfK\n\n",\
    particles, species, r_1, r_2, val_1, val_2, rp_d, diff_c, box_len, dt,\
    max_eq_time, max_t_steps, msd_steps, gr_steps,flag, eps_r, temp);
}

bool read_file_atom_pos(){
    
    FILE* param = fopen("input_mono_rcp.dat", "r");

    if(param == (FILE*)NULL){
        return false;
    }
    else{
        
    }
}