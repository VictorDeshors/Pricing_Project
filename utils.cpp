#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>

// Standard normal cumulative distribution function (CDF)
double normalCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
}

// Black-Scholes function to price European call or put options
double blackScholesPrice(double S, double K, double T, double r, double sigma, std::string optionType) {
    // Calculate d1 and d2
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    if (optionType == "Call") {
        // Call option price: S * N(d1) - K * exp(-r * T) * N(d2)
        return S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
    }
    else if (optionType == "Put") {
        // Put option price: K * exp(-r * T) * N(-d2) - S * N(-d1)
        return K * std::exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
    }
    else{
        throw std::runtime_error("Unknown option type.");

    }
}

// Discretization fct
std::vector<double> getGrid(double init, double step, int length){
    std::vector<double> time_grid(length, 0.0);
    for (int i = 0; i < length; i++) {
        time_grid[i] = init + i * step;
    }
    return  time_grid;
}

void printVector(std::vector<double> vec){
    for (double elem: vec){
        std::cout << elem << std::endl;
    }
}
