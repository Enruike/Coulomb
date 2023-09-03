#include "functions.hpp"

double box_muller(double num1, double num2){

    double result;

    result = cos(2.0 * M_PI * num2) * sqrt(-2 * log(num1));

    return result;
}