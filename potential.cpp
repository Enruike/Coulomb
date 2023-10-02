#include"potential.hpp"
#include"constants.hpp"
#include"functions.hpp"

int main(){

    FILE* electric_energy;
    FILE* repulsive_energy;
    FILE* mean_square_displacement;
    FILE* electrolyte_movie;
    FILE * rhor_file;
    FILE * gr_file;

    if(!read_parameters()){
        printf("No parameters file!\n");
        exit(1);
    }
   
    Bjerrum_len = pow(q_electron, 2) / (4. * M_PI * eps_r * epsilon_zero * K_boltzmann * temp);
    Bjerrum_red = Bjerrum_len / 1e-10;
    diff_c_red = diff_c / 1.e-10;
    diff_c_red2 = diff_c_red / 1.e-10;

    //calculating grid dimension. It starts from 0.
    dim_gr = rint(.5 * box_len / delta_gr);

    printf("Bjerrum length is %.6e\n", Bjerrum_len);
    printf("Reduced Bjerrum length is %.6e\n", Bjerrum_red);
    printf("Grid dimension %d\n", dim_gr);

    /* VERY IMPORTANT 
        Macroflag will be set to false till macroion implementation
        is ready. */
    macro_flag = false;

    //Number of atoms per specie
    int * atoms_per_specie = (int*)malloc(species * sizeof(int));
    
    if(macro_flag){
        /* nothing till now */
    }
    else{
        for(int i = 0; i < species; i++){
            atoms_per_specie[i] = num_particles / species;
        }
    }

    /* Condiciones de electroneutralidad */
    double charge_val;
    if(species == 2){

        charge_val = val_1 * atoms_per_specie[0] + val_2 * atoms_per_specie[1];
        
        if(charge_val != 0.){
            printf("Electronegativity condition is not met!\n");
            exit(1);
        }
    }

    char rhor_file_name[100];
    char gr_file_name[100];
    double XR[dim_gr] = { 0. };
    double XRP[dim_gr] = { 0. };
    double bin_vol[dim_gr] = { 0. };

    for(int i = 0; i < dim_gr; i++){
        XR[i] = (i + 0.5) * delta_gr;
        XRP[i] = XR[i] - XR[0];
    }

    for(int i = 0; i < dim_gr; i++){
        if(i != dim_gr - 1){
            bin_vol[i] = 4. / 3. * M_PI * (pow(XRP[i + 1], 3) - pow(XRP[i], 3));
        }
        else{
            bin_vol[i] = 4. / 3. * M_PI * (pow((dim_gr + 0.5) * delta_gr - XR[0], 3));
        }
    }

    //double HR[species][species][dim_gr] = { 0. };
    //double GR[species][species][dim_gr] = { 0. };
    //double RHOR[species][species][dim_gr] = { 0. };
    
    double ***HR = new double**[species];
    double ***HR_temp = new double ** [species];
    double ***RHOR = new double **[species];
    double ***GR = new double **[species];

    for(int i = 0; i < species; i++){

        HR[i] = new double * [species];
        HR_temp[i] = new double * [species];
        RHOR[i] = new double * [species];
        GR[i] = new double *[species];

        for(int j = 0; j < species; j++){

            HR[i][j] = new double[dim_gr];
            HR_temp[i][j] = new double[dim_gr];
            RHOR[i][j] = new double[dim_gr];
            GR[i][j] = new double[dim_gr];

            for(int k = 0; k < dim_gr; k++){
                HR[i][j][k] = 0.;
                HR_temp[i][j][k] = 0.;
                RHOR[i][j][k] = 0.;
                GR[i][j][k] = 0.;
            }

        }
    }

    /*
    printf("value for HR[1][0][50] %lf\n", HR[1][0][50]);
    HR[1][0][50] = 3.56;
    printf("value for HR[1][0][50] %lf\n", HR[1][0][50]);
    */

    //printf("rint(3.5) %lf, int(3.5) %d\n", rint(3.5), (int)3.5);

    //printf("K boltzman is %.1e\n", K_boltzmann);
    //printf("Pi is %lf\n\n", M_PI);

    //Particle *particles = new Particle[num_particles];
    
    if(!read_file_atom_pos()){
        printf("No positions file found!\n");
        exit(1);
    }

    printf("%lf\t%lf\t%lf\n", positions[0 * 3 + 0], positions[0 * 3 + 1], positions[0 * 3 + 2]);
    printf("%lf\t%lf\t%lf\n", positions[99 * 3 + 0], positions[99 * 3 + 1], positions[99 * 3 + 2]);

    very_initial_positions = (double*)malloc(num_particles * 3 * sizeof(double));
    /* Positions x, y, and z are stored in array named positions. 
        We are gonna multiply these positions by reduced length box at once. */
    for(int i = 0; i < num_particles; i++){

        positions[i * 3 + 0] = (positions[i * 3 + 0] - 0.5) * box_len;
        positions[i * 3 + 1] = (positions[i * 3 + 1] - 0.5) * box_len;
        positions[i * 3 + 2] = (positions[i * 3 + 2] - 0.5) * box_len;

        very_initial_positions[i * 3 + 0] = positions[i * 3 + 0];
        very_initial_positions[i * 3 + 1] = positions[i * 3 + 1];
        very_initial_positions[i * 3 + 2] = positions[i * 3 + 2];
    
        /*
        particles[i].pos_x = positions[i * 3 + 0];
        particles[i].pos_y = positions[i * 3 + 1];
        particles[i].pos_z = positions[i * 3 + 2];
        */

    }

    printf("Diameter %lf\n", diameter);
    printf("Max diameter %lf\n", diameter * box_len);

    printf("%lf\t%lf\t%lf\n", positions[0 * 3 + 0], positions[0 * 3 + 1], positions[0 * 3 + 2]);
    printf("%lf\t%lf\t%lf\n", positions[(num_particles - 1) * 3 + 0], positions[99 * 3 + 1], positions[99 * 3 + 2]);
    

    vol_frac = num_particles * (M_PI / 6.) * pow(diameter, 3); //      !!! Check the use of volume fraction

    printf("Volume fraction: %lf\n", vol_frac);

    char_array = (char*)malloc(num_particles * sizeof(char));
    radius_array = (double*)malloc(num_particles * sizeof(double));
    valence_array = (double*)malloc(num_particles * sizeof(double));
    species_array = (short int*)malloc(num_particles * sizeof(short int));

    for(int i = 0; i < num_particles; i++){
        if(i < num_particles / species){
            //particles[i].specie = 1;
            //particles[i].radius = r_1;
            radius_array[i] = r_1;
            //particles[i].valence = val_1;
            valence_array[i] = val_1;
            char_array[i] = 'A';
            species_array[i] = 0;
        }
        else{
            //particles[i].specie = 2;
            //particles[i].radius = r_2;
            radius_array[i] = r_2;
            //particles[i].valence = val_2;
            valence_array[i] = val_2;
            char_array[i] = 'B';
            species_array[i] = 1;
        }
        //particles[i].infinite = 0;
    }

    /* Random variables */

    unsigned seed = 101013;
    double random_number;
    double ran_num1, ran_num2;

    //seed the random number generator
    std::default_random_engine generator(seed);

    //Create a uniform real distribution
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    printf("\n\t##### random number generator #####\n");
    for(int i = 0; i < 3; i++){
        random_number = distribution(generator);
        printf("%lf\n", random_number);
    }

    printf("\t##### Box Muller generator #####\n");
    
    for(int i = 0; i < 10; i++){
        ran_num1 = distribution(generator);
        ran_num2 = distribution(generator);
        printf("%lf\n", box_muller(ran_num1, ran_num2));
    }

    printf("\n\t\t#####\n");


    /* Generador de archivo movie.xml */
    //HOOMD_xml_generator();

    //
    double val_rc = 0.;
    double val_el = 0.;
    double e_rc = 0.;
    double e_el = 0.;
    double dxij, dyij, dzij, rij;
    double x_pos, y_pos, z_pos;
    
    double xi, yi, zi;
    
    for(int indx = 0; indx < num_particles; indx++){
        
        xi = positions[indx * 3 + 0];
        yi = positions[indx * 3 + 1];
        zi = positions[indx * 3 + 2];        

        for(int i = 0; i < num_particles; i++){

            if(i != indx){

                periodic_distance(xi, yi, zi, x_pos, y_pos, z_pos, &positions[i * 3]);

                dxij = xi - x_pos;
                dyij = yi - y_pos;
                dzij = zi - z_pos;
                rij = sqrt(pow(dxij, 2) + pow(dyij, 2) + pow(dzij, 2));
                e_rc += erc(rij, radius_array[indx], radius_array[i]);
                e_el += eew(rij, valence_array[indx], valence_array[i]);
                //printf("i : %d\n", i);
                //printf("e_el %.6e\n", 0.5 * e_el / num_particles);
                
            }
        }
        
        /* Acumuladores para la energía total que sienten las partículas */
        val_rc += 0.5 * (e_rc / num_particles);
        val_el += 0.5 * (e_el / num_particles);
        //printf("e_el %lf \n", 0.5 * e_el / num_particles);

        /* Debido a que no uso una función para calcular las distancias
           del i-ésimo átomo con respecto del j-ésimo, debo reiniciar el
           acumulador para las energías del repulsive core y la potencial.*/
        e_rc = 0; //Aquí reseteo la energía del repulsive core.
        e_el = 0; //Aquí pongo a cero el potencial eléctrico.

        //printf("val_el %.6e\n", val_el);
    }
    

    printf("Repulsive Core Energy %.6e\n", val_rc);
    printf("Electrostatic Energy %.6e\n\t\t#####\n\n", val_el);

    /* Creando archivos de energía eléctrica y repulsiva */
    electric_energy = fopen("electric_energy.out", "w");
    repulsive_energy = fopen("repulsive_energy.out", "w");

    /* First column corresponds to iteration times time step 
       de second column is the energy value. */
    fprintf(repulsive_energy, "0\t\t%.6e\n", val_rc);
    fprintf(electric_energy, "0\t\t%.6e\n", val_el);

    fclose(repulsive_energy);
    fclose(electric_energy);

    /* Creating MSD file */
    mean_square_displacement = fopen("mean_square_displacement.out", "w");
    fclose(mean_square_displacement);

    /* Electrolyte movie file */
    electrolyte_movie = fopen("electrolyte_movie.xyz", "w");
    fclose(electrolyte_movie);

    /* ##### This is the most important part #### */

    double * new_positions = (double*)malloc(num_particles * 3 * sizeof(double));
    cells = (short int*)malloc(num_particles * 3 * sizeof(short int));
    double msd_val = 0.;
    double tau = 1.;

    //Inicializando valores

    for(int i = 0; i < num_particles; i++){
        for(int j = 0; j < 3; j++){
            cells[i * 3 + j] = 0;
            new_positions[i * 3 + j] = 0.;
        }
    }

    //printf("%d, %d, %d\n", cells[0], cells[1], cells[2]);
    //printf("%d, %d, %d\n", cells[99 * 3], cells[99 * 3 + 1], cells[99 * 3 + 2]);
    //printf("dt is %.e\n", dt);
    //max_t_steps

    printf("##### Equilibration initialized #####\n");
    printf("\t# Set to %d steps #\n", min_eq_steps);
    for(int iter = 0; iter < max_time_steps; iter++){
        
        for(int indx = 0; indx < num_particles; indx++){
            
            if(std::isnan(new_pos_function(indx, dt, &positions[indx * 3], &new_positions[indx * 3], &cells[indx * 3]))){
                printf("Moving particles function error!\n");
                printf("iter step %d indx %d\n", iter, indx);
                exit(1);
            }
        }

        for(int i = 0; i < num_particles; i++){

            //if((i == 0 || i == num_particles - 1) && (iter + 1) % 5000 == 0) printf("%lf\t%lf\t%lf\n",positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
            //if(i == 99 && ) printf("old %lf\t%lf\t%lf\n",positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
            positions[i * 3 + 0] = new_positions[i * 3 + 0];
            positions[i * 3 + 1] = new_positions[i * 3 + 1];
            positions[i * 3 + 2] = new_positions[i * 3 + 2];
            //if(i == 99) printf("new %lf\t%lf\t%lf\n",positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
            //if((i == 0 || i == num_particles - 1) && (iter + 1) % 5000 == 0) printf("%lf\t%lf\t%lf\n",positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);
        
        }

        //if((iter + 1) % (msd_steps) == 0 && iter != 0){
        if((iter + 1) % msd_steps == 0){    
    
            /* Mean Square Displacement implementation */
            for(int i = 0; i < num_particles; i++){
                msd_val += pow(positions[i * 3 + 0] + (double)cells[i * 3 + 0] * box_len - very_initial_positions[i * 3 + 0], 2)\
                + pow(positions[i * 3 + 1] + (double)cells[i * 3 + 1] * box_len - very_initial_positions[i * 3 + 1], 2)\
                + pow(positions[i * 3 + 2] + (double)cells[i * 3 + 2] * box_len - very_initial_positions[i * 3 + 2], 2);
            }

            msd_val /= num_particles;

            mean_square_displacement = fopen("mean_square_displacement.out", "a");
            
            fprintf(mean_square_displacement, "%.6e\t\t%.6e\n", (iter + 1) * dt, msd_val);

            fclose(mean_square_displacement);

            msd_val = 0.;

            /* Esta sección fue copiada directamente de la parte inicial.
               Aquí se calcula la energía total electrostática y de núcleo
               repulsivo. */
            
            val_rc = 0.;
            val_el = 0.;

            for(int indx = 0; indx < num_particles; indx++){
            //for(int indx = 0; indx < 1; indx++){

                xi = positions[indx * 3 + 0];
                yi = positions[indx * 3 + 1];
                zi = positions[indx * 3 + 2];        

                for(int i = 0; i < num_particles; i++){

                    if(i != indx){

                        periodic_distance(xi, yi, zi, x_pos, y_pos, z_pos, &positions[i * 3]);

                        dxij = xi - x_pos;
                        dyij = yi - y_pos;
                        dzij = zi - z_pos;
                        rij = sqrt(pow(dxij, 2) + pow(dyij, 2) + pow(dzij, 2));
                        if(std::isnan(e_rc += erc(rij, radius_array[indx], radius_array[i]))){
                            printf("An error has ocurred!\n");
                            printf("iter step %d indx %d\n", iter, indx);
                            exit(1);
                        }
                        e_el += eew(rij, valence_array[indx], valence_array[i]);

                        //printf("e_rc %lf, e_el %lf\n",0.5 * (e_rc / num_particles), 0.5 * (e_el / num_particles));
                        
                    }
                    
                }
                //printf("e_el temp %lf\n", 0.5 * (e_el / num_particles));
                val_rc += e_rc * 0.5 / num_particles;
                val_el += e_el * 0.5 / num_particles;
               
                //printf("e_rc %lf, e_el %lf\n",0.5 * (e_rc / num_particles), 0.5 * (e_el / num_particles));

                e_rc = 0;
                e_el = 0;

            }
            
            //printf("#### Iteration %d ####\n", iter + 1);
            //printf("val_rc %.6e\n", val_rc);
            //printf("val_el %.6e\n", val_el);

            //printf("%d, %d, %d\n", cells[0], cells[1], cells[2]);
            //printf("%d, %d, %d\n", cells[99 * 3], cells[99 * 3 + 1], cells[99 * 3 + 2]);

            electric_energy = fopen("electric_energy.out", "a");
            repulsive_energy = fopen("repulsive_energy.out", "a");

            fprintf(repulsive_energy, "%.6e\t\t%.6e\n", (iter + 1) * dt, val_rc);
            fprintf(electric_energy, "%.6e\t\t%.6e\n", (iter + 1) * dt, val_el);

            fclose(repulsive_energy);
            fclose(electric_energy);

        }

        if((iter + 1) % histo_steps == 0 && (iter + 1) > min_eq_steps){
            
            if(tau == 1){
                printf("... Equilibration is done!\n\n");
            }

            /* Movie */
            
            electrolyte_movie = fopen("electrolyte_movie.xyz", "a");
            fprintf(electrolyte_movie,"%d\n", num_particles);
            fprintf(electrolyte_movie, "Electrolyte\n");
            for(int i = 0; i < num_particles; i++){
                fprintf(electrolyte_movie, "%c\t%.6e\t%.6e\t%.6e\n", char_array[i], positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);   
            }
            fclose(electrolyte_movie);

            /* Histograma */

            histogram_hr_tau(num_particles, positions, species_array, HR_temp);
            
            for(int i = 0; i < species; i++){
                for(int j = 0; j < species; j++){
                    for(int k = 0; k < dim_gr; k++){
                        HR[i][j][k] += HR_temp[i][j][k];
                        HR_temp[i][j][k] = 0.;
                    }
                }
            }

            if((int)tau > 0 && (int)tau % tau_steps == 0){

                printf("##### %d\u03c4 iteration #####\n", (int)tau);

                calculate_rhor_gr(RHOR, GR, HR, tau, atoms_per_specie, bin_vol);

                snprintf(rhor_file_name, sizeof(rhor_file_name), "%d_rhor.out", (int)tau);
                snprintf(gr_file_name, sizeof(gr_file_name), "%d_gr.out", (int)tau);
                
                rhor_file = fopen(rhor_file_name, "w");
                gr_file = fopen(gr_file_name, "w");

                fclose(rhor_file);
                fclose(gr_file);

                write_gr_rhor(gr_file, gr_file_name, (int)tau, species, XR, GR);
                write_gr_rhor(rhor_file, rhor_file_name, (int)tau, species, XR, RHOR);

                /* Valores a 0 */
                for(int i = 0; i < species; i++){
                    for(int j = 0; j < species; j++){
                        for(int k = 0; k < dim_gr; k++){
                            RHOR[i][j][k] = 0.;
                            GR[i][j][k] = 0.;
                        }
                    }
                }
            }

            tau += 1.;
        }
    }

    /*for(int i = 0; i < species; i++){
        for(int j = 0; j < species; j++){
            for(int k = 0; k < dim_gr; k++){
                printf("%lf ", HR[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }*/

    /* Free dynamic memory arrays */
    free(positions);
    free(new_positions);
    free(very_initial_positions);
    free(species_array);
    free(cells);
    free(char_array);
    free(radius_array);
    free(valence_array);
    free(atoms_per_specie);

    //free HR_temp triple-pointer
    for(int i = 0; i < species; i++){
        for(int j = 0; j < species; j++){
            delete[] HR[i][j];
            delete[] HR_temp[i][j];
            delete[] RHOR[i][j];
            delete[] GR[i][j];
        }
        delete[] HR_temp[i];
        delete[] HR[i];
        delete[] RHOR[i];
        delete[] GR[i];
    }
    delete[] HR_temp;
    delete[] HR;
    delete[] RHOR;
    delete[] GR;

    /* Destroy class dynamic memory array */
    //delete[] particles;

    return 0;
    //printf("x for particle 3 is %lf\n", particle[2].pos_x); //Segmentation fault because we erase the array. No more memory.
}