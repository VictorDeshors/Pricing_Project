#include <iostream>
#include <algorithm>
#include <format>
#include <math.h>
#include <ctime>
#include "pde_pricer.h"
#include "simulations.h"
#include "utils.h"
#include <fstream>
#include <vector>
#include <string>



int main()
{
    /////// PARAMS ///////
    // Timing variables
    clock_t start, end;
    int size_test = 100000;
    double T = 0.5; // time to expiry*/
    double atm_vol = 0.16; // ATM volatility
    double K = 45.0; // strike of the option
    double S0 = 50.0; // spot price on the day of pricing
    double r = 0.03, q = 0.0; // risk-free rate, assumed constant and repo rate including div yield
    
    // Schema param
    double theta = 0.5; // choosing Crank-Nicholson scheme

    // Time grid params
    int n = 50; // nb of steps to discretize time grid
    std::vector<double> time_grid = getGrid(0, T / n, n + 1);

    // Space grid param
    int m = 100; // nb of steps to discretize the space grid
    double lambda = 4.0;
    double s_m = std::max(S0 * (1 - lambda * sqrt(T) * atm_vol), 0.0);
    double s_M = S0 * (1 + lambda * sqrt(T) * atm_vol);
    std::vector<double> space_grid = getGrid(s_m, (s_M - s_m) / m, m + 1);

    /////// CONDITIONS ///////
    // Boundary conditions
    std::vector<std::vector<double>> bound_cond(n + 1,
                                                std::vector<double>(2, 0.0)); // should be of size (n rows, 2 cols)
    for (int i = 0; i < n + 1; i++) {
        bound_cond[i][0] = std::max(s_m - K * exp(-(r - q) * time_grid[i]), 0.0);
        bound_cond[i][1] = s_M - K * exp(-(r - q) * (T - time_grid[i]));
    }

    // Terminal conditions
    std::vector<double> term_cond(m + 1, 0.0);  // should be of size (m rows, 1 col)
    for (int i = 0; i < m + 1; i++) {
        term_cond[i] = std::max(space_grid[i] - K, 0.0);
    }

    // PDE params
    Matrix A(std::vector<std::vector<double>>(n + 1, std::vector<double>(m + 1, -r)));
    Matrix B(std::vector<std::vector<double>>(n + 1, std::vector<double>(m + 1, r - q)));
    Matrix C(std::vector<std::vector<double>>(n + 1, std::vector<double>(m + 1, atm_vol * atm_vol / 2)));
    Matrix D(std::vector<std::vector<double>>(n + 1, std::vector<double>(m + 1, 0.0)));
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            double space_elem = space_grid[j];
            B.data[i][j] *= space_elem;
            C.data[i][j] *= space_elem * space_elem;
        }
    }

    /////// SOLVING ///////
    // Solving with BS
    start = clock();
    double bs_price = blackScholesPrice(S0, K, T, r, atm_vol, "Call");
    end = clock();
    std::cout << "B.S. price of: $" << std::format("{:.2f}", bs_price) << " found in "
              << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;


    // Solving with the PDE
    start = clock();
    PDEPricer pde_pricer = PDEPricer(T, n, m, space_grid, bound_cond, term_cond, A, B, C, D, theta);
    double pde_price = pde_pricer.getOptionPrice(S0);
    end = clock();
    std::cout << "Finite Difference price of: $" << std::format("{:.2f}", pde_price) << " found in "
                << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;





    



    /////// PLOTTING
    std::cout << "\n Plotting prices..." << std::endl;
    std::vector<std::pair<double, double>> results;
    
    // Effect of varying time to maturity on the price of the option with PDE
    for (double T_prime = 0.1; T_prime <= 2.0; T_prime += 0.1) {
        // Update time grid
        std::vector<double> time_grid = getGrid(0, T_prime / n, n + 1);

        // Update Boundary conditions
        std::vector<std::vector<double>> bound_cond(n + 1, std::vector<double>(2, 0.0)); // should be of size (n rows, 2 cols)
        for (int i = 0; i < n + 1; i++) {
            bound_cond[i][0] = std::max(s_m - K * exp(-(r - q) * time_grid[i]), 0.0);
            bound_cond[i][1] = s_M - K * exp(-(r - q) * (T_prime - time_grid[i]));
        }

        // Solving with the PDE
        start = clock();
        PDEPricer pde_pricer = PDEPricer(T_prime, n, m, space_grid, bound_cond, term_cond, A, B, C, D, theta);
        double pde_price = pde_pricer.getOptionPrice(S0);
        end = clock();
        std::cout << "Finite Difference price of: $" << std::format("{:.2f}", pde_price) << " for T_prime = " << T_prime << " found in "
                << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;

        // Store the result
        results.emplace_back(T_prime, pde_price);
    }
    // Save results for time to maturity to a CSV file
    std::ofstream file_time_to_maturity("../prices_to_plot/time_to_maturity.csv");
    file_time_to_maturity << "TimeToMaturity,OptionPrice\n";
    for (const auto& result : results) {
        file_time_to_maturity << result.first << "," << result.second << "\n";
    }
    file_time_to_maturity.close();
    


    //// Effect of varying volatility on the price of the option with PDE
    results.clear();
    // Loop through different values of volatility
    for (double atm_vol_prime = 0.1; atm_vol_prime <= 0.5; atm_vol_prime += 0.01) {
        // Update C matrix
        C = Matrix(std::vector<std::vector<double>>(n + 1, std::vector<double>(m + 1, atm_vol_prime * atm_vol_prime / 2)));
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                double space_elem = space_grid[j];
                C.data[i][j] *= space_elem * space_elem;
            }
        }
        // Solving with the PDE
        start = clock();
        PDEPricer pde_pricer = PDEPricer(T, n, m, space_grid, bound_cond, term_cond, A, B, C, D, theta);
        double pde_price = pde_pricer.getOptionPrice(S0);
        end = clock();
        std::cout << "Finite Difference price of: $" << std::format("{:.2f}", pde_price) << " for ATM Volatility = " << atm_vol_prime << " found in "
                << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;

        // Store the result
        results.emplace_back(atm_vol_prime, pde_price);
    }
    // Save results for volatility to a CSV file
    std::ofstream file_volatility("../prices_to_plot/volatility.csv");
    file_volatility << "ATMVolatility,OptionPrice\n";
    for (const auto& result : results) {
        file_volatility << result.first << "," << result.second << "\n";
    }
    file_volatility.close();




    // Effect of varying time steps on the price of the option with PDE
    results.clear();
    // Loop through different values of time steps
    for (int time_step = 50; time_step <= 150; time_step += 5) {
        std::cout << "Testing with time_step = " << time_step << std::endl;

        // Update time grid
        std::vector<double> time_grid = getGrid(0, T / time_step, time_step + 1);

        // Update PDE params
        Matrix A_prime(std::vector<std::vector<double>>(time_step + 1, std::vector<double>(m + 1, -r)));
        Matrix B_prime(std::vector<std::vector<double>>(time_step + 1, std::vector<double>(m + 1, r - q)));
        Matrix C_prime(std::vector<std::vector<double>>(time_step + 1, std::vector<double>(m + 1, atm_vol * atm_vol / 2)));
        for (int i = 0; i < time_step + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                double space_elem = space_grid[j];
                B_prime.data[i][j] *= space_elem;
                C_prime.data[i][j] *= space_elem * space_elem;
            }
        }
        Matrix D_prime(std::vector<std::vector<double>>(time_step + 1, std::vector<double>(m + 1, 0.0)));

        // Update Boundary conditions
        std::vector<std::vector<double>> bound_cond_prime(time_step + 1,
                                                        std::vector<double>(2, 0.0)); 
        for (int i = 0; i < time_step + 1; i++) {
            bound_cond_prime[i][0] = std::max(s_m - K * exp(-(r - q) * time_grid[i]), 0.0);
            bound_cond_prime[i][1] = s_M - K * exp(-(r - q) * (T - time_grid[i]));
        }

        // Solving with the PDE
        start = clock();
        PDEPricer pde_pricer = PDEPricer(T, time_step, m, space_grid, bound_cond_prime, term_cond, A_prime, B_prime, C_prime, D_prime, theta);
        double pde_price = pde_pricer.getOptionPrice(S0);
        end = clock();
        std::cout << "Finite Difference price of: $" << std::format("{:.2f}", pde_price) << " for n = " << time_step << " found in "
                << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;

        // Store the result
        results.emplace_back(time_step, pde_price);
    }
    // Compute the theoretical price using Black-Scholes formula
    bs_price = blackScholesPrice(S0, K, T, r, atm_vol, "Call");
    // Save results for time steps to a CSV file
    std::ofstream file_time_steps("../prices_to_plot/time_steps.csv");
    file_time_steps << "TimeSteps,OptionPrice,BS_price\n";
    for (const auto& result : results) {
        file_time_steps << result.first << "," << result.second << "," << bs_price << "\n";
    }
    file_time_steps.close();







    /// Effect of varying space steps on the price of the option with PDE
    results.clear();
    // Loop through different values of space steps
    for (int space_step = 90; space_step <= 300; space_step += 30) {
        std::cout << "Testing with space_step = " << space_step << std::endl;

        // Update space grid
        std::vector<double> space_grid = getGrid(s_m, (s_M - s_m) / space_step, space_step + 1);

        // Update Boundary conditions
        std::vector<std::vector<double>> bound_cond(n + 1, std::vector<double>(2, 0.0)); // should be of size (n rows, 2 cols)
        for (int i = 0; i < n + 1; i++) {
            bound_cond[i][0] = std::max(s_m - K * exp(-(r - q) * time_grid[i]), 0.0);
            bound_cond[i][1] = s_M - K * exp(-(r - q) * (T - time_grid[i]));
        }

        // Update Terminal conditions
        std::vector<double> term_cond(space_step + 1, 0.0);  // should be of size (space_step rows, 1 col)
        for (int i = 0; i < space_step + 1; i++) {
            term_cond[i] = std::max(space_grid[i] - K, 0.0);
        }

        // Update PDE params
        Matrix A(std::vector<std::vector<double>>(n + 1, std::vector<double>(space_step + 1, -r)));
        Matrix B(std::vector<std::vector<double>>(n + 1, std::vector<double>(space_step + 1, r - q)));
        Matrix C(std::vector<std::vector<double>>(n + 1, std::vector<double>(space_step + 1, atm_vol * atm_vol / 2)));
        Matrix D(std::vector<std::vector<double>>(n + 1, std::vector<double>(space_step + 1, 0.0)));

        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < space_step + 1; j++) {
                double space_elem = space_grid[j];
                B.data[i][j] *= space_elem;
                C.data[i][j] *= space_elem * space_elem;
            }
        }

        // Solving with the PDE
        start = clock();
        PDEPricer pde_pricer = PDEPricer(T, n, space_step, space_grid, bound_cond, term_cond, A, B, C, D, theta);
        double pde_price = pde_pricer.getOptionPrice(S0);
        end = clock();
        std::cout << "Finite Difference price of: $" << std::format("{:.2f}", pde_price) << " for m = " << space_step << " found in "
                << (end - start) * 1000.0 / CLOCKS_PER_SEC << "ms" << std::endl;

        // Store the result
        results.emplace_back(space_step, pde_price);
    }
    // Save results for space steps to a CSV file
    // Compute the theoretical price using Black-Scholes formula
    bs_price = blackScholesPrice(S0, K, T, r, atm_vol, "Call");
    std::ofstream file_space_steps("../prices_to_plot/space_steps.csv");
    file_space_steps << "SpaceSteps,OptionPrice,BS_price\n";
    for (const auto& result : results) {
        file_space_steps << result.first << "," << result.second << "," << bs_price << "\n";
    }
    file_space_steps.close();







    #include <fstream>
    /// PART USED FOR MC AND QMC METHODS
    std::cout << "MC and QMC Prices: " << std::endl;
    std::vector<double> params = getGrid(100, 100, size_test / 100);
    int nb_sim_mc = 500;
    bs_price = blackScholesPrice(S0, K, T, r, atm_vol, "Call");

    std::vector<int> nb_sims(size_test / 100, 0);
    std::vector<double> vec_mc_prices(size_test / 100, 0.0);
    std::vector<double> vec_qmc_prices(size_test / 100, 0.0);
    for (int i = 0; i < size_test / 100; i++) {
        int nb_sim = params[i];
        nb_sims[i] = nb_sim;

        double mc_price = MCEuCallPricer(S0, K, T, r, q, atm_vol, nb_sim);
        vec_mc_prices[i] = mc_price;

        QMCEuCallPricer qmc_pricer = QMCEuCallPricer(S0, K, T, r, q, atm_vol, nb_sim);
        double qmc_price = qmc_pricer.get_qmc_price();
        vec_qmc_prices[i] = qmc_price;
        //    std::cout << "Done for nb sim: " << nb_sim << std::endl; 
    }


    // Save results to a CSV file
    std::ofstream file("../prices_to_plot/mc_qmc_prices.csv");
    file << "NbSim,MCPrice,QMCPrice,BS_price\n";
    for (size_t i = 0; i < vec_mc_prices.size(); ++i) {
        file << nb_sims[i] << "," << vec_mc_prices[i] << "," << vec_qmc_prices[i] << "," << bs_price << "\n";
    }
    file.close();

    return 1;
}
