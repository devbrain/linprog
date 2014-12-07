#ifndef __PROBLEM_DESCRIPTION_HPP__
#define __PROBLEM_DESCRIPTION_HPP__

#include <vector>

struct problem_description
{
    unsigned int m; // number of constraints
    unsigned int n; // number of vars

    std::vector <unsigned> basic; // indices of basic variables |basic| = m
    std::vector <unsigned> non_basic; // indices of basic variable |non_basic| = n

    std::vector <double> b; // constraint boundaries values |b| = m
    std::vector <double>   a; // m x n matrix of constraints

    double z0; // objective function value
    std::vector <double> c; // coefficients

    problem_description(const char* fname); // load from file
};


#endif
