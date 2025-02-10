#include <matplot/matplot.h>
#include <eulerCauchy.h>
#include <improvedEulerCauchy.h>
#include <rungeKutta4.h>
#include <algorithm>
#include <cmath>


std::vector<double> computeErrorNorm(std::vector<double> y_exact, std::vector<double> y) {

    std::vector<double> error;

    for (size_t i = 0; i < y.size(); ++i) {
        error.push_back(std::abs(y_exact[i] - y[i]));
    }

    return error;
}


int main(){
    
    using namespace matplot;

    // Setting up ODE
    std::function<double(double, double)> odeFunction = [](double x, double y){
        return x*x*x - y;
    };
    std::pair<double, double> xInit = {0, 0};
    double xEnd = 2;
    double h = 0.1;

    // Exact Solution
    std::vector<double> xExact = iota(xInit.first, h, xEnd);
    std::vector<double> yExact(xExact.size());
    std::transform(xExact.begin(), xExact.end(), yExact.begin(), [](double x) {
        return pow(x, 3) - 3 * pow(x, 2) + 6 * x - 6 + 6 * exp(-x);
    });

    // ODE Integration
    std::vector<std::pair<double, double>> eulerCauchySol = eulerCauchySolver(odeFunction, xInit, xEnd, h);
    std::vector<std::pair<double, double>> improvedEulerCauchySol = improvedEulerCauchySolver(odeFunction, xInit, xEnd, h);
    std::vector<std::pair<double, double>> rungeKutta4Sol = rungeKutta4Solver(odeFunction, xInit, xEnd, h);

    // Extract Solution
    std::vector<double> xEulerCauchy(eulerCauchySol.size()), yEulerCauchy(eulerCauchySol.size());
    std::vector<double> xImprovedEulerCauchy(improvedEulerCauchySol.size()), yImprovedEulerCauchy(improvedEulerCauchySol.size());
    std::vector<double> xRungeKutta4(rungeKutta4Sol.size()), yRungeKutta4(rungeKutta4Sol.size());

    std::transform(eulerCauchySol.begin(), eulerCauchySol.end(), xEulerCauchy.begin(),
                   [](const std::pair<double, double>& p) { return p.first; });
    std::transform(eulerCauchySol.begin(), eulerCauchySol.end(), yEulerCauchy.begin(),
                   [](const std::pair<double, double>& p) { return p.second; });
    std::transform(improvedEulerCauchySol.begin(), improvedEulerCauchySol.end(), xImprovedEulerCauchy.begin(),
                   [](const std::pair<double, double>& p) { return p.first; });
    std::transform(improvedEulerCauchySol.begin(), improvedEulerCauchySol.end(), yImprovedEulerCauchy.begin(),
                   [](const std::pair<double, double>& p) { return p.second; });
    std::transform(rungeKutta4Sol.begin(), rungeKutta4Sol.end(), xRungeKutta4.begin(),
                   [](const std::pair<double, double>& p) { return p.first; });
    std::transform(rungeKutta4Sol.begin(), rungeKutta4Sol.end(), yRungeKutta4.begin(),
                   [](const std::pair<double, double>& p) { return p.second; });

    // Plotting
    figure();
    gcf()->size(800, 600);
    
    hold(on);
    plot(xExact, yExact);
    plot(xEulerCauchy, yEulerCauchy);
    plot(xImprovedEulerCauchy, yImprovedEulerCauchy);
    plot(xRungeKutta4, yRungeKutta4);
    hold(off);

    xlabel("x");
    ylabel("u");
    legend({"Exact", "Euler-Cauchy", "Improved Euler-Cauchy", "Runge-Kutta 4"});

    save("results/odePlot.jpg");

    figure();
    gcf()->size(800, 600);

    hold(on);
    loglog(xExact, computeErrorNorm(yExact, yEulerCauchy));
    loglog(xExact, computeErrorNorm(yExact, yImprovedEulerCauchy));
    loglog(xExact, computeErrorNorm(yExact, yRungeKutta4));
    hold(off);

    xlabel("x");
    ylabel("Error norm");
    legend({"Euler-Cauchy", "Improved Euler-Cauchy", "Runge-Kutta 4"});
    

    save("results/errorPlot.jpg");
    
    return 0;
}