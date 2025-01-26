#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include "simulations.h"

// Classic Monte Carlo method- a function suffices here
double MCEuCallPricer(double S0, double K, double T, double r, double q, double sigma, int numSimulations) {
    // Initialize random number generator for uniform distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0); // Standard normal distribution

    double totalPayoff = 0.0;
    // Simulating numSimulations final paysoffs
    for (int i = 0; i < numSimulations; i++) {
        double Z = d(gen);
        double ST = S0 * std::exp((r - q - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
        double payoff = std::max(ST - K, 0.0);
        totalPayoff += payoff;
    }
    // Discount the avg payoff to present value
    double optionPrice = std::exp(-r * T) * (totalPayoff / numSimulations);
    return optionPrice;
}

// Quasi-MC method; wrapping all the components within a class for ease of use
QMCEuCallPricer::QMCEuCallPricer(double S0, double K, double T, double r, double q, double sigma, int numSimulations):
S0(S0), K(K), T(T), r(r), q(q), sigma(sigma), numSimulations(numSimulations){}

void QMCEuCallPricer::initialize_sobol_sequence() {
    // Initialize direction numbers for the first dimension
    std::vector<unsigned int> directionNumbers(32, 0);
    for (int i = 0; i < 32; i++) {
        directionNumbers[i] = 1u << (31 - i);  // bit-shifting instead of std::pow
    }
    // Precompute normalization factor
    const double norm_factor = 1.0 / (1ULL << 32);  // 1 / 2^32
    std::vector<double> sobolSq(numSimulations, 0.0);
    unsigned int y = 0;  // Initial value y_0 = 0
    for (int i = 0; i < numSimulations; i++) {
        // Find the rightmost zero bit in i (i + 1 is used since sequence is 1-based)
        int j = 0;
        int temp = i + 1;
        while (temp & 1) {
            temp >>= 1;
            j++;
        }
        // XOR with the corresponding direction number
        y ^= directionNumbers[j];
        // Normalization
        sobolSq[i] = static_cast<double>(y) * norm_factor;
    }
    shuffleVector(sobolSq);
    this->sobolSequence = sobolSq;
}


void QMCEuCallPricer::initialize_quasi_random_normal_variables(){
    // Having initialized our quasi random sequence, we use the Box Muller algorithm to get a vector of normal variables
    std::vector<double> normal_sq(numSimulations, 0.0);
    int rest_eucli_num_sim = (numSimulations % 2);
    // If uneven number of simulations, I reuse a random elem of my Sobol sequence to compute the last elem
    for (int j=0; j < numSimulations - 1 - rest_eucli_num_sim; j=j+2){
        double u1 = sobolSequence[j]; double u2 = sobolSequence[j+1];
        double z1 = sqrt(-2 * log(u1)) * cos(2* M_PI * u2);
        double z2 = sqrt(-2 * log(u1)) * sin(2* M_PI * u2);
        normal_sq[j] = z1; normal_sq[j + 1] = z2;
    }
    if (rest_eucli_num_sim == 1){
        std::random_device rd;
        std::mt19937 gen(rd());  // Mersenne Twister PRNG
        std::uniform_int_distribution<int> dist(0, numSimulations - 1);
        int random_index = dist(gen);
        double u1 = sobolSequence[numSimulations - 1]; double u2 = sobolSequence[random_index];
        double z1 = sqrt(-2 * log(u1)) * cos(2* M_PI * u2);
        normal_sq[numSimulations - 1] = z1;
    }
    this->normalSequence = normal_sq;
}

double QMCEuCallPricer::get_qmc_price(){
    initialize_sobol_sequence();
    initialize_quasi_random_normal_variables();
    if (numSimulations < 1){
        throw std::runtime_error("Please enter a higher number of simulations.");

    }
    double totalPayoff = 0.0;
    // Simulating numSimulations final paysoffs
    for (int i = 0; i < numSimulations; i++) {
        double Z = normalSequence[i];
        double ST = S0 * std::exp((r - q - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
        double payoff = std::max(ST - K, 0.0);
        totalPayoff += payoff;
    }
    // Discount the avg payoff to present value
    double optionPrice = std::exp(-r * T) * (totalPayoff / numSimulations);
    return optionPrice;
}

void QMCEuCallPricer::shuffleVector(std::vector<double>& vec) {
    std::random_device rd;   // Random seed generator
    std::mt19937 gen(rd());  // Mersenne Twister PRNG
    std::shuffle(vec.begin(), vec.end(), gen);
}