#include "functions.hpp"

double box_muller(double num1, double num2){

    double result;

    result = cos(2.0 * M_PI * num2) * sqrt(-2 * log(num1));

    return result;
}

void HOOMD_xml_generator(){

    FILE *movie;
    movie = fopen("movie.xml", "w");

    fprintf(movie, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(movie, "<hoomd_xml version=\"1.7\">\n");
    fprintf(movie, "<configuration time_step=\"1\" dimensions=\"3\" natoms=\"%d\">\n", num_particules);
    fprintf(movie, "<box lx=\"30.0\" ly=\"30.0\" lz=\"30.0\" xy=\"0\" xz=\"0\" yz=\"0\"/>\n");
    fprintf(movie, "<position num=\"%d\">\n", num_particules);
    
    for(int i = 0; i < num_particules; i++){

        fprintf(movie, "%lf\t%lf\t%lf\n", positions[i * 3 + 0], positions[i * 3 + 1], positions[i * 3 + 2]);

    }

    fprintf(movie, "</position>\n");
    fprintf(movie, "<mass num=\"0\">\n");
    fprintf(movie, "</mass>\n");
    fprintf(movie, "<charge num=\"0\">\n");
    fprintf(movie, "</charge>\n");
    fprintf(movie, "<diameter num=\"0\">\n");
    fprintf(movie, "</diameter>\n");
    fprintf(movie, "<type num=\"0\">\n");

    for(int i = 0; i < num_particules; i++){
        fprintf(movie, "%c\n", char_array[i]);
    }

    fprintf(movie, "</type>\n");
    fprintf(movie, "<body num=\"0\">\n");
    fprintf(movie, "</body>\n");
    fprintf(movie, "<bond num=\"0\">\n");
    fprintf(movie, "</bond>\n");
    fprintf(movie, "<angle num=\"0\">\n");
    fprintf(movie, "</angle>\n");
    fprintf(movie, "<dihedral num=\"0\">\n");
    fprintf(movie, "</dihedral>\n");
    fprintf(movie, "<improper num=\"0\">\n");
    fprintf(movie, "</improper>\n");
    fprintf(movie, "</configuration>\n");
    fprintf(movie, "</hoomd_xml>\n");

    fclose(movie);

}

double energy_rc_i_all(int indx, int np_total){

    double e_rc = 0.;
    double dxij, dyij, dzij, rij;
    
    double x_array[np_total], y_array[np_total], z_array[np_total];
    double x_current, y_current, z_current;

    double xi = positions[indx * 3 + 0];
    double yi = positions[indx * 3 + 1];
    double zi = positions[indx * 3 + 2];
    
    double x_diff = 0., y_diff = 0., z_diff = 0.;

    for(int i = 0; i < np_total; i++){

        x_current = positions[i * 3 + 0];
        y_current = positions[i * 3 + 1];
        z_current = positions[i * 3 + 2];

        x_diff = x_current - xi;
        y_diff = y_current - yi;
        z_diff = z_current - zi;

        //Loop for X values
        if(x_diff > half_box){
            x_array[i] = x_current - box_len;
        }
        else if(x_diff < - half_box){
            x_array[i] = x_current + box_len;
        }
        else{
            x_array[i] = x_current;
        }

        //Loop for Y values
        if(y_diff > half_box){
            y_array[i] = y_current - box_len;
        }
        else if(y_diff < - half_box){
            y_array[i] = y_current + box_len;
        }
        else{
            y_array[i] = y_current;
        }

        //Loop for Z values
        if(z_diff > half_box){
            z_array[i] = z_current - box_len;
        }
        else if(z_diff < - half_box){
            z_array[i] = z_current + box_len;
        }
        else{
            z_array[i] = z_current;
        }
    }

    for(int i = 0; i < np_total; i++){
        if(i != indx){

            dxij = xi - x_array[i];
            dyij = yi - y_array[i];
            dzij = zi - z_array[i];
            rij = sqrt(pow(dxij, 2) + pow(dyij, 2) + pow(dzij, 2));
            e_rc += erc(rij, radius_array[indx], radius_array[i]);

        }
    }

    return 0.5 * (e_rc / np_total);
}

double erc(double r, double ri, double rj){

    double dij;
    double x;
    double d_rc = 1.0; //Sigma for the rc potential.
    double res;

    dij = ri + rj - d_rc;

    if(r <= dij){
        printf("r : %lf\n dij : %lf\n", r, dij);
        printf("Infinite Electrostatic Energy!\n");
        exit(1);        
    }
    else{
        x = r - dij;
    }

    if(r < (dij + pow(2. * d_rc, 1. / 6.))){
        res = 4. * (pow(d_rc / x, 12) - pow(d_rc / x, 6)) + 1.;
    }
    else{
        res = 0.;
    }

    return res;

}

double energy_el_i_all(int indx, int np_total){

    double e_rc = 0.;
    double dxij, dyij, dzij, rij;
    
    double x_array[np_total], y_array[np_total], z_array[np_total];
    double x_current, y_current, z_current;

    double xi = positions[indx * 3 + 0];
    double yi = positions[indx * 3 + 1];
    double zi = positions[indx * 3 + 2];
    
    double x_diff = 0., y_diff = 0., z_diff = 0.;

    for(int i = 0; i < np_total; i++){

        x_current = positions[i * 3 + 0];
        y_current = positions[i * 3 + 1];
        z_current = positions[i * 3 + 2];

        x_diff = x_current - xi;
        y_diff = y_current - yi;
        z_diff = z_current - zi;

        //Loop for X values
        if(x_diff > half_box){
            x_array[i] = x_current - box_len;
        }
        else if(x_diff < - half_box){
            x_array[i] = x_current + box_len;
        }
        else{
            x_array[i] = x_current;
        }

        //Loop for Y values
        if(y_diff > half_box){
            y_array[i] = y_current - box_len;
        }
        else if(y_diff < - half_box){
            y_array[i] = y_current + box_len;
        }
        else{
            y_array[i] = y_current;
        }

        //Loop for Z values
        if(z_diff > half_box){
            z_array[i] = z_current - box_len;
        }
        else if(z_diff < - half_box){
            z_array[i] = z_current + box_len;
        }
        else{
            z_array[i] = z_current;
        }
    }

    for(int i = 0; i < np_total; i++){
        if(i != indx){

            dxij = xi - x_array[i];
            dyij = yi - y_array[i];
            dzij = zi - z_array[i];
            rij = sqrt(pow(dxij, 2) + pow(dyij, 2) + pow(dzij, 2));
            e_rc += eew(rij, valence_array[indx], valence_array[i]);

        }
    }

    //printf("%lf\t", 0.5 * (e_rc / np_total));
    return 0.5 * (e_rc / np_total);

}

double eew(double r, double vi, double vj){
    return Bjerrum_red * ( (vi * vj) / r );
}