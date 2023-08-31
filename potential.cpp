#include"potential.hpp"
#include"constants.hpp"

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

    Particule *particules = new Particule[num_particules];
    
    if(!read_file_atom_pos()){
        printf("No positions file found!\n");
        exit(1);
    }

    /* Positions x, y, and z are stored in array named positions. 
        We are gonna multiply these positions by reduced length box at once. */
    for(int i = 0; i < num_particules; i++){
        particules[i].pos_x = positions[i * 3 + 0] * box_len;
        particules[i].pos_y = positions[i * 3 + 1] * box_len;
        particules[i].pos_z = positions[i * 3 + 2] * box_len;
    }

    printf("Diameter %lf\n", diameter);
    printf("Max diameter %lf\n", diameter * box_len);

    vol_frac = num_particules * (M_PI / 6.) * pow(diameter, 3); //      !!! Check the use of volume fraction

    printf("Volume fraction: %lf\n", vol_frac);

    for(int i = 0; i < num_particules; i++){
        if(i < num_particules / species){
            particules[i].specie = 1;
            particules[i].radius = r_1;
            particules[i].valence = val_1;
        }
        else{
            particules[i].specie = 2;
            particules[i].radius = r_2;
            particules[i].valence = val_2;
        }
        particules[i].infinite = 0;
    }

    free(positions);
    delete[] particules;
    return 0;
    //printf("x for particle 3 is %lf\n", particle[2].pos_x); //Segmentation fault because we erase the array. No more memory.
}