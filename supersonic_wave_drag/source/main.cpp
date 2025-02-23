#include <solver.h>
#include <matplot/matplot.h>
#include <chrono>

int main(){
    using namespace matplot;

    // Initializing vectors
    std::vector<double> M_vec{1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    std::vector<double> cd_num_vec;
    std::vector<double> cd_exact_vec;

    // Airfoil constant
    const double e_m = 0.1;

    
    // Computing Drag Coefficient
    for (double M : M_vec){
        cd_num_vec.push_back(full_domain(M));
        cd_exact_vec.push_back(double (16.0/3.0 * e_m * e_m / (std::sqrt(M*M - 1))));
    }

    // Start timing
    auto start_full_domain = std::chrono::high_resolution_clock::now();

    full_domain(M_vec[0]);

    // Finish timing
    auto end_full_domain = std::chrono::high_resolution_clock::now();

    // Start timing
    auto start_sliding_window = std::chrono::high_resolution_clock::now();

    sliding_window(M_vec[0]);

    // Finish timing
    auto end_sliding_window = std::chrono::high_resolution_clock::now();


    // Outputting
    std::chrono::duration<double> elapsed_full_domain = end_full_domain - start_full_domain;
    std::chrono::duration<double> elapsed_sliding_window = end_sliding_window - start_sliding_window;
    std::cout << "Elapsed time (Full Domain): " << elapsed_full_domain.count() << " seconds\n";
    std::cout << "Elapsed time (Sliding Window): " << elapsed_sliding_window.count() << " seconds\n";
    for (int i = 0; i < static_cast<int>(M_vec.size()); ++i){
        std::cout << "@ M = " << M_vec[i] << " , cd = " << cd_num_vec[i] << std::endl;
    }

    // Plotting
    figure();
    gcf()->size(800, 600);
    
    hold(on);
    plot(M_vec, cd_num_vec);
    plot(M_vec, cd_exact_vec);
    hold(off);
    
    xlabel("M");
    ylabel("cd");
    legend({"Numerical", "Exact"});

    save("results/cd_plot.jpg");

    return 0;
}