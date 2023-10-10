#include "read_parameters.hpp"

bool read_parameters(){

    file_pos_gen = 0;
    
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
        fscanf(param, "Temperature = %lf #Temperature in Kelvin(K)\n\n", &temp);
        fscanf(param, "### Macroion parameters ###\n\n");
        fscanf(param, "Number of macroions = %d #0: no macroion implementation\n", &macro_num);

        if(macro_num != 0){
            fscanf(param, "Valence = %lf\n", &macro_valence);
            fscanf(param, "Radius = %lf\n", &macro_radius);
            fscanf(param, "Repulsive core distance = %lf\n", &macro_rc);
            fscanf(param, "Position file generator = %d #0 if you have the positions file\n", &file_pos_gen);
        }
    }
    
    char* flag;

    if(disol == 0){
        flag = "False";
    }
    else{
        flag = "True";
    }

    if(species != 2 ){
        printf("Code just works with 2 species!\n");
        printf("*** BAD TERMINATION ***\n");
        exit(1);
    }

    printf("**** PARAMETERS FILE ****\n\n\
    Number of particles = %d\n\
    Number of species = %d\n\
    Radius of species 1 = %.2lf\u00c5\n\
    Radius of species 2 = %.2lf\u00c5\n\
    Valence of species 1 = %+.2lf\n\
    Valence of species 2 = %+.2lf\n\
    Repulsive core distance = %.2lf\u00c5\n\
    Diffusion coefficient = %.2e\n\
    Repulsive core sigma = %.2lf\n\
    Length of the box = %.2lf\u00c5\n\
    Delta grid = %.2lf\n\
    Delta time = %.2es\n\
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

    if(macro_num != 0) printf("**** MACROION PARAMETERS ****\n\n\
    Number of macroions = %d\n\
    Valence = +%.2lf\n\
    Radius = %.2lf\u00c5\n\
    Repulsive Core Distance = %.2lf\u00c5\n\
    Position file generator = %d\n\n",\
    macro_num, macro_valence, macro_radius, macro_rc, file_pos_gen);

    fclose(param);

    return true;
    
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

    fclose(param);

    return true;
}

//std::default_random_engine generator(101013);
//std::default_random_engine generator(static_cast<unsigned int>(std::time(nullptr)));
//std::uniform_real_distribution<double> distribution(0., 1.);

bool pos_macro_gen(){

    double pos[3] = { 0. };
    double macro_1[3] = { 0. };
    double macro_2[3] = { 0. };
    double dis = 0., dis_m1 = 0., dis_m2 = 0.;
    double xij = 0., yij = 0., zij = 0.;
    double m1xij = 0., m1yij = 0., m1zij = 0.;
    double m2xij = 0., m2yij = 0., m2zij = 0.;
    bool overlap = true;
    double radius = 0.; //radius for the bigger species;

    /* Como nota, cabe mencionar que no tomaremos en cuenta la distancia del repulsive core para esta implementacion. */
    if(r_1 > r_2){
        radius = r_1;
    }
    else{
        radius = r_2;
    }

    positions = (double*)malloc(num_particles * 3 * sizeof(double));
    
    positions[(num_particles - 2) * 3 + 0] = macro_1[0] = .5;
    positions[(num_particles - 2) * 3 + 1] = macro_1[1] = 0.5;
    positions[(num_particles - 2) * 3 + 2] = macro_1[2] = 0.5;
    positions[(num_particles - 1) * 3 + 0] = macro_2[0] = 0.75;
    positions[(num_particles - 1) * 3 + 1] = macro_2[1] = 0.75;
    positions[(num_particles - 1) * 3 + 2] = macro_2[2] = 0.75;

    /* for(int i = 0; i < 3; i++){
        macro_1[i] = (macro_1[i] - 0.5) * box_len;
        macro_2[i] = (macro_2[i] - 0.5) * box_len;
    } */

    for(int i = (num_particles - macro_num); i < num_particles; i++){
        printf("x: %.10lf y:%.10lf z:%.10lf\n", positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
    }
    printf("m1x: %lf m1y:%lf m1z:%lf\n", macro_1[0], macro_1[1], macro_1[2]);
    printf("m2x: %lf m2y:%lf m2z:%lf\n", macro_2[0], macro_2[1], macro_2[2]);

    printf("num particles is %d\n", num_particles);
    //printf("macro num %d\n", macro_num);

    for(int i = 0; i < (num_particles - macro_num); i++){
        
        printf("iter i:%d\n", i);
        pos_gen(pos);
        //printf("x: %.10lf y:%.10lf z:%.10lf\n", pos[0], pos[1], pos[2]);
        while (overlap){
            for(int j = 0; j <= i; j++){     
            
                if(j != i){
                    printf("i:%d j:%d ", i, j);
                    
                        
                    xij = (positions[j * 3 + 0] - pos[0]) * box_len;
                    yij = (positions[j * 3 + 1] - pos[1]) * box_len;
                    zij = (positions[j * 3 + 2] - pos[2]) * box_len;

                    m1xij = (macro_1[0] - pos[0]) * box_len;
                    m1yij = (macro_1[1] - pos[1]) * box_len;
                    m1zij = (macro_1[2] - pos[2]) * box_len;
                    
                    m2xij = (macro_2[0] - pos[0]) * box_len;
                    m2yij = (macro_2[1] - pos[1]) * box_len;
                    m2zij = (macro_2[2] - pos[2]) * box_len;

                /*  //length conversion
                    xij = xij * box_len;
                    yij = yij * box_len;
                    zij = zij * box_len;

                    m1xij = m1xij * box_len;
                    m1yij = m1yij * box_len;
                    m1zij = m1zij * box_len;
                    
                    m2xij = m2xij * box_len;
                    m2yij = m2yij * box_len;
                    m2zij = m2zij * box_len; */


                    //distance ij
                    dis = sqrt(pow(xij, 2) +pow(yij, 2) + pow(zij, 2));
                    dis_m1 = sqrt(pow(m1xij, 2) +pow(m1yij, 2) + pow(m1zij, 2));
                    dis_m2 = sqrt(pow(m2xij, 2) +pow(m2yij, 2) + pow(m2zij, 2));

                    printf("dis:%lf dis_m1:%lf dis_m2:%lf\n", dis, dis_m1, dis_m2);

                    if(dis > (2. * radius) && dis_m1 > (macro_radius + radius) && dis_m2 > (macro_radius + radius)){
                        overlap = false;
                    }
                    else{
                        printf("Overlapping ocurred! New position generated.\n");
                        pos_gen(pos);
                        j = -1;
                    }

                }
                       
            }

            if(i == 0){
                overlap = false;
            }
            
        }
        
        positions[i * 3 + 0] = pos[0];
        positions[i * 3 + 1] = pos[1];
        positions[i * 3 + 2] = pos[2];

        overlap = true; 
        
    }

    printf("positions array\n");
    for(int i = 0; i < num_particles; i++){
        printf("x: %.10lf y:%.10lf z:%.10lf\n", positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
    }

}

bool pos_macro_read(){

    FILE * macro_pos = fopen("pos_macro.in", "r");

    if(macro_pos == (FILE*)NULL){
        printf("No pos_macro.in file found!\n");
        return false;
    }
    else{

    }

    fclose(macro_pos);

    return true;

}