#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <cmath>

// Takes in Mach number and outputs drag coefficient
double sliding_window(double M);
double full_domain(double M);

#endif