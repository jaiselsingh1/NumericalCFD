#ifndef IMPROVEDEULERCAUCHY_H
#define IMPROVEDEULERCAUCHY_H

#include <vector>
#include <functional>

std::vector<std::pair<double, double>> improvedEulerCauchySolver(
    std::function<double(double, double)> f, // Input function
    std::pair<double, double> xInit, // Initial conditions
    double xEnd, // Solve until x = xEnd
    double h // Step size
);

#endif