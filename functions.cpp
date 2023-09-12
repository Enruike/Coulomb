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
    fprintf(movie, "<configuration time_step=\"1\" dimensions=\"3\" natoms=\"%d\">\n", num_particles);
    fprintf(movie, "<box lx=\"30.0\" ly=\"30.0\" lz=\"30.0\" xy=\"0\" xz=\"0\" yz=\"0\"/>\n");
    fprintf(movie, "<position num=\"%d\">\n", num_particles);
    
    for(int i = 0; i < num_particles; i++){

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

    for(int i = 0; i < num_particles; i++){
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
    double res;

    dij = ri + rj - d_rc;

    if(r <= dij){
        printf("r : %lf\n dij : %lf\n", r, dij);
        printf("Infinite Electrostatic Energy!\n");
        return std::nan("");        
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

void periodic_distance(double xi, double yi, double zi,
    double& x_pos, double& y_pos, double& z_pos, double* positions){

        double x_current, y_current, z_current;
        double x_diff, y_diff, z_diff;

        x_current = positions[0];
        y_current = positions[1];
        z_current = positions[2];

        x_diff = x_current - xi;
        y_diff = y_current - yi;
        z_diff = z_current - zi;

        //Loop for X values
        if(x_diff > half_box){
            x_pos = x_current - box_len;
        }
        else if(x_diff < -half_box){
            x_pos = x_current + box_len;
        }
        else{
            x_pos = x_current;
        }

        //Loop for Y values
        if(y_diff > half_box){
            y_pos = y_current - box_len;
        }
        else if(y_diff < -half_box){
            y_pos = y_current + box_len;
        }
        else{
            y_pos = y_current;
        }

        //Loop for Z values
        if(z_diff > half_box){
            z_pos = z_current - box_len;
        }
        else if(z_diff < -half_box){
            z_pos = z_current + box_len;
        }
        else{
            z_pos = z_current;
        }

    }

double new_pos_function(int indx, double dt, double * indx_positions, double new_pos[3], int * cells){

    double xi, yi, zi;
    double fx, fy, fz;
    double gx, gy, gz;
    double x_j, y_j, z_j;
    double fx_rc_KBT = 0., fy_rc_KBT = 0., fz_rc_KBT = 0.;
    double fx_ew_KBT = 0., fy_ew_KBT = 0., fz_ew_KBT = 0.;
    double frc_r = 0.;
    double few_r = 0.;
    double dxij = 0., dyij = 0., dzij = 0.;
    double rij;

    xi = indx_positions[0];
    yi = indx_positions[1];
    zi = indx_positions[2];
    //printf("index %d\n ", indx);

    for(int i = 0; i < num_particles; i++){

        
        if(i != indx){
            
            periodic_distance(xi, yi, zi, x_j, y_j, z_j, &positions[i * 3]);
            //printf("x:%lf y:%lf z:%lf\n", new_pos[0], new_pos[1], new_pos[2]);
            //printf("xi %lf yi %lf zi %lf \n", xi, yi, zi);
            //printf("x_pos %lf y_pos %lf z_pos %lf\n", x_pos, y_pos, z_pos);
            dxij = xi - x_j;
            dyij = yi - y_j;
            dzij = zi - z_j;

            rij = sqrt(pow(dxij, 2) + pow(dyij, 2) + pow(dzij, 2));
            //printf("rij %lf, dxij %lf, dyij %lf, dzij %lf\n", rij, dxij, dyij, dzij);
            frc_r = frc(rij, radius_array[indx], radius_array[i]);
            //printf("repulsive %.6e\n", frc_r);
            fx_rc_KBT += frc_r * dxij / rij;
            fy_rc_KBT += frc_r * dyij / rij;
            fz_rc_KBT += frc_r * dzij / rij;
            
            few_r = few0(rij, valence_array[indx], valence_array[i]);
            //printf("electric %.6e\n", few_r);
            fx_ew_KBT += few_r * dxij / rij;
            fy_ew_KBT += few_r * dyij / rij;
            fz_ew_KBT += few_r * dzij / rij;

        }
    }

    //printf("repulsive sum %.6e\n", frc_r);
    //printf("electric sum %.6e\n", few_r);

    fx = fx_rc_KBT + fx_ew_KBT;
    fy = fy_rc_KBT + fy_ew_KBT;
    fz = fz_rc_KBT + fz_ew_KBT;

    //printf("fx %.6e fy %.6e fz %.6e\n", fx, fy, fz);

    gx = random_muller();
    gy = random_muller();
    gz = random_muller();

    //printf("gx %lf gy %lf gz %lf\n", gx, gy, gz);

    new_pos[0] = xi + fx * diff_c_red * dt + gx * sqrt(2. * diff_c_red2 * dt);
    new_pos[1] = yi + fy * diff_c_red * dt + gy * sqrt(2. * diff_c_red2 * dt);
    new_pos[2] = zi + fz * diff_c_red * dt + gz * sqrt(2. * diff_c_red2 * dt);
    //printf("xi is %.6e\n", xi);
    //printf("Part 1 is %.6e\n", (diff_c / 1.e-10) * dt);
    //printf("diff_C 2 is %lf\n", diff_c_red2);
    //printf("Part 2 is %.6e\n", gx * sqrt(2. * diff_c_red2 * dt));
    //printf("new pos x should be %.6e\n", xi + fx * (diff_c / 1.e-10) * dt + gx * sqrt(2. * dt * diff_c / 1.e-20));
    //printf("pos inside function %lf\t%lf\t%lf\n", new_pos[0], new_pos[1], new_pos[2]);
    

    if(std::isnan(mi_after_move(indx, new_pos[0], new_pos[1], new_pos[2], &cells[indx * 3]))){
        printf("\n!!!\n\tProblem after moving particles!\n");
        printf("\t\tdt is too large\n!!!\t\n\n");
        printf("original positions x:%lf y:%lf z:%lf\n", xi, yi, zi);
        printf("x:%lf y:%lf z:%lf\n", new_pos[0], new_pos[1], new_pos[2]);
        return std::nan("");
    }

    //printf("after mi move x:%lf y:%lf z:%lf\n", new_pos[0], new_pos[1], new_pos[2]);

    return 0.;

}

double frc(double rij, double ri, double rj){

    double dij;
    double x;
    double s;
    double res;

    dij = ri + rj - d_rc;

    if(rij <= dij){
        s = 1.000000001 * dij - dij;
    }
    else{
        s = rij - dij;
    }

    double term = dij + d_rc * pow(2., 1./6.);

    if(rij < term){
        res = 1. / (1.e-10 * d_rc) * (48. * pow(d_rc / s, 13.) - 24. * pow(d_rc / s, 7.));
    }
    else{
        res = 0.;
    }

    return res;
}

double few0(double rij, double vi, double vj){
    
    double res;
    //printf("Bjerrum %lf\n", Bjerrum_len);
    //printf("vi %lf vj %lf\n", vi, vj);
    //printf("rij %lf\n", rij);
    //printf("pow %lf\n", pow(rij, 2));
    res = Bjerrum_len * (vi * vj) / pow(rij, 2);
    //printf("res %.6e\n", res);
    return res;

}

std::default_random_engine generator(101013);
//std::default_random_engine generator(static_cast<unsigned int>(std::time(nullptr)));
std::uniform_real_distribution<double> distribution(0., 1.);

double random_muller(){

    double res;
    double num1, num2;

    num1 = distribution(generator);
    num2 = distribution(generator);

    //printf("num1 %lf num2 %lf\n", num1, num2);

    res = cos(2.0 * M_PI * num2) * sqrt(-2 * log(num1));

    return res;

}

double mi_after_move(int indx, double& x, double& y, double& z, int * cell){

    //printf("indx %d begining mi_after_function x:%lf y:%lf z:%lf\n",indx, x, y, z);

    if(x > half_box && x <= 3. * half_box){
        x -= box_len;
        cell[0] += 1;
    }
    else if(x < - half_box && x >= -3. * half_box){
        x += box_len;
        cell[0] -= 1;
    }
    else if(x > 3. * half_box || x < -3. * half_box){
        printf("\nx problem\n");
        printf("x:%lf\n", x);
        return std::nan("");  
    }

    if(y > half_box && y <= 3. * half_box){
        y -= box_len;
        cell[1] += 1;
    }
    else if(y < - half_box && y >= -3. * half_box){
        y += box_len;
        cell[1] -= 1;
    }
    else if(y > 3. * half_box || y < -3. * half_box){
        printf("\ny problem\n");
        printf("y:%lf\n", y);
        return std::nan("");  
    }

    if(z > half_box && z <= 3. * half_box){
        z -= box_len;
        cell[2] += 1;
    }
    else if(z < -half_box && z >= -3. * half_box){
        z += box_len;
        cell[2] -= 1;
    }
    else if(z > 3. * half_box || z < -3. * half_box){
        printf("z problem\n");
        printf("z:%lf\n", z);
        return std::nan("");  
    }

    //printf("indx %d end mi_after_function x:%lf y:%lf z:%lf\n", indx, x, y, z);

    return 0.;
}

void calculate_energies(){
    
}