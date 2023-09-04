#include"potential.hpp"
#include"constants.hpp"
#include"functions.hpp"

int main(){
    if(!read_parameters()){
        printf("No parameters file!\n");
        exit(1);
    }
   
    Bjerrum_len = pow(q_electron, 2) / (4. * M_PI * eps_r * K_boltzmann * temp);
    Bjerrum_red = Bjerrum_len / 1e-10;

    printf("Bjerrum length is %.6e\n", Bjerrum_len);
    printf("Reduced Bjerrum length is %.6e\n", Bjerrum_red);
    //printf("K boltzman is %.1e\n", K_boltzmann);
    //printf("Pi is %lf\n\n", M_PI);

    //Particule *particules = new Particule[num_particules];
    
    if(!read_file_atom_pos()){
        printf("No positions file found!\n");
        exit(1);
    }

    printf("%lf\t%lf\t%lf\n", positions[0 * 3 + 0], positions[0 * 3 + 1], positions[0 * 3 + 2]);
    printf("%lf\t%lf\t%lf\n", positions[99 * 3 + 0], positions[99 * 3 + 1], positions[99 * 3 + 2]);

    /* Positions x, y, and z are stored in array named positions. 
        We are gonna multiply these positions by reduced length box at once. */
    for(int i = 0; i < num_particules; i++){

        positions[i * 3 + 0] = (positions[i * 3 + 0] - 0.5) * box_len;
        positions[i * 3 + 1] = (positions[i * 3 + 1] - 0.5) * box_len;
        positions[i * 3 + 2] = (positions[i * 3 + 2] - 0.5) * box_len;

        /*
        particules[i].pos_x = positions[i * 3 + 0];
        particules[i].pos_y = positions[i * 3 + 1];
        particules[i].pos_z = positions[i * 3 + 2];
        */

    }

    printf("Diameter %lf\n", diameter);
    printf("Max diameter %lf\n", diameter * box_len);

    vol_frac = num_particules * (M_PI / 6.) * pow(diameter, 3); //      !!! Check the use of volume fraction

    printf("Volume fraction: %lf\n", vol_frac);

    char_array = (char*)malloc(num_particules * sizeof(char));
    radius_array = (double*)malloc(num_particules * sizeof(double));
    valence_array = (double*)malloc(num_particules * sizeof(double));

    for(int i = 0; i < num_particules; i++){
        if(i < num_particules / species){
            //particules[i].specie = 1;
            //particules[i].radius = r_1;
            radius_array[i] = r_1;
            //particules[i].valence = val_1;
            valence_array[i] = val_1;
            char_array[i] = 'A';
        }
        else{
            //particules[i].specie = 2;
            //particules[i].radius = r_2;
            radius_array[i] = r_2;
            //particules[i].valence = val_2;
            valence_array[i] = val_2;
            char_array[i] = 'B';
        }
        //particules[i].infinite = 0;
    }

    /* Random variables */

    unsigned seed = 101013;
    double random_number;
    double ran_num1, ran_num2;

    //seed the random number generator
    std::default_random_engine generator(seed);

    //Create a uniform real distribution
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    printf("\t##### random number generator #####\n");
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

    printf("\t\t#####\n");

    HOOMD_xml_generator();

    //
    double val_rc = 0.;
    for(int indx = 0; indx < num_particules; indx++){
        val_rc += energy_rc_i_all(indx, num_particules);
    }

    printf("val_rc %.6e\n", val_rc);

    double val_el = 0.;

    for(int indx = 0; indx < num_particules; indx++){
        val_el += energy_el_i_all(indx, num_particules);
    }
    
    printf("val_el %.6e\n", val_el);

    /* Free dynamic memory arrays */
    free(positions);
    free(char_array);
    free(radius_array);
    free(valence_array);

    /* Destroy class dynamic memory array */
    //delete[] particules;

    return 0;
    //printf("x for particle 3 is %lf\n", particle[2].pos_x); //Segmentation fault because we erase the array. No more memory.
}