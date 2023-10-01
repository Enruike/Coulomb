#include "read_parameters.hpp"

bool read_parameters(){

    FILE* param = fopen("param.in", "r");

	if (param == (FILE*)NULL) {
		return false;
	}
    else{
        fscanf(param, "**** PARAMETERS FILE ****\n\n");
        fscanf(param, "Number of particles = %d\n", &num_particles);
        fscanf(param, "Number of species = %d\n", &species);
        fscanf(param, "Radius of species 1 = %lf\n", &r_1);
        fscanf(param, "Radius of species 2 = %lf\n", &r_2);
        fscanf(param, "Valence of species 1 = %lf\n", &val_1);
        fscanf(param, "Valence of species 2 = %lf\n", &val_2);
        fscanf(param, "Repulsive core distance = %lf\n", &rp_d);
        fscanf(param, "Diffusion coefficient = %lf\n", &diff_c);
        fscanf(param, "Repulsive core sigma = %lf #Hardness\n", &d_rc);
        fscanf(param, "Length of the box = %lf #Length in Angstroms\n", &box_len);
        fscanf(param, "Delta grid = %lf #1 space over 5 splits = 0.2\n", &delta_gr);
        fscanf(param, "Delta time = %lf #Time step in seconds\n", &dt);
        fscanf(param, "Min equilibration time steps = %d #Min time for equilibration\n", &min_eq_steps);
        fscanf(param, "Max time steps = %d #Max number of steps\n", &max_time_steps);
        fscanf(param, "Energy steps = %d #Saving steps for Energy and MSD file\n", &msd_steps);
        fscanf(param, "Histogram steps = %d #Steps for calculating histogram\n", &histo_steps);
        fscanf(param, "tau steps = %d #Steps for calculating g(r) and rho(r)\n0", &tau_steps);
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
    Diffusion coefficient = %.2e\n\
    Repulsive core sigma = %.2lf\n\
    Length of the box = %.2lf\u00c5\n\
    Delta grid = %.2lf\n\
    Delta time = %.2e s\n\
    Min equilibration time steps = %d\n\
    Max time steps = %d\n\
    Energy & MSD saving steps = %d\n\
    Histogram steps = %d\n\
    \u03c4 steps = %d\n\
    Infinite disolution = %s\n\
    Epsilon r = %.2lf\n\
    Temperature = %.2lfK\n\n",\
    num_particles, species, r_1, r_2, val_1, val_2, rp_d, diff_c, d_rc,box_len, delta_gr, dt,\
    min_eq_steps, max_time_steps, msd_steps, histo_steps, tau_steps,flag, eps_r, temp);

    half_box = box_len / 2.;
}

bool read_file_atom_pos(){
    
    FILE* param = fopen("input_mono_rcp.dat", "r");

    if(param == (FILE*)NULL){
        return false;
    }
    else{

        positions = (double*)malloc(num_particles * 3 * sizeof(double));

        for(int i = 0; i < 6; i++){
            
            if(i == 3){
                fscanf(param, "%lf\n", &diameter);
            }
            else{
                fscanf(param, "%*[^\n]\n");
            }
            
        }

        for(int i = 0; i < num_particles; i++){

            fscanf(param, "%lf %lf %lf\n", &positions[i * 3 + 0], &positions[i * 3 + 1], &positions[i * 3 + 2]);
            
        }
       /*  for(int i = 0; i < num_particles; i++){
            printf("x: %lf \t y: %lf \t z: %lf\n", positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
        } */
    }
}