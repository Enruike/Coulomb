#include"potential.hpp"
#include"constants.hpp"

int main(){
    if(!read_parameters()){
        printf("No parameters file!\n");
        exit(1);
    }
    else{
        Bjerrum_len = pow(q_electron, 2) / (4. * M_PI * eps_r * K_boltzmann * temp);
        Bjerrum_red = Bjerrum_len / 1e-10;
        printf("Working fine till now!\n");
        printf("Bjerrum length is %.6e\n", Bjerrum_len);
        printf("Reduced Bjerrum length is %.6e\n", Bjerrum_red);
        printf("K boltzman is %.1e\n", K_boltzmann);
        printf("Pi is %lf\n", M_PI);
    }
}