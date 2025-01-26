#include "pde_pricer.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

PDEPricer::PDEPricer(
    double T, int n, int m, std::vector<double> space_grid,
    const std::vector<std::vector<double>>& bound_cond, const std::vector<double>& term_cond,
    const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& D, double theta
):T(T), n(n), m(m), space_grid(space_grid), bound_cond(bound_cond), term_cond(term_cond),
A(A), B(B), C(C), D(D), theta(theta){}

Matrix PDEPricer::getPi(int timestep_i, double delta_x){
    // Getting the params
    double delta_t = T / n;
    std::vector<double> a_i = A.data[timestep_i];
    std::vector<double> b_i = B.data[timestep_i];
    std::vector<double> c_i = C.data[timestep_i];
    // Filling the Pi matrix
    std::vector<std::vector<double>> P_i(m - 1, std::vector<double>( m - 1, 0.0));
    for (int j = 0; j < m - 1; j++) {
        double current_a = a_i[j + 1]; double current_b = b_i[j + 1]; double current_c = c_i[j + 1];
        // Main diagonal
        P_i[j][j] = current_a - (1.0 / delta_t + 2.0 * theta * current_c / (delta_x * delta_x));
        // Lower diagonal
        if (j > 0) {
            P_i[j][j - 1] = -current_b / (2.0 * delta_x) + theta * current_c / (delta_x * delta_x);
        }
        // Upper diagonal
        if (j < m - 2) {
            P_i[j][j + 1] = current_b / (2.0 * delta_x) + theta * current_c / (delta_x * delta_x);
        }
    }
    return Matrix(P_i);
}

Matrix PDEPricer::getQi(int timestep_i, double delta_x){
    // Getting the params
    double delta_t = T / n;
    std::vector<double> c_i = C.data[timestep_i];

    // Filling the Qi matrix
    std::vector<std::vector<double>> Q_i(m - 1, std::vector<double>(m - 1, 0.0));
    for (int j = 0; j < m - 1; j++) {
        double current_c =  c_i[j + 1];
        // Main diagonal
        Q_i[j][j] = (1.0 / delta_t) - 2.0 * (1.0 - theta) * current_c / (delta_x * delta_x);
        // Lower diagonal
        if (j > 0) {
            Q_i[j][j - 1] = (1.0 - theta) * current_c / (delta_x * delta_x);
        }
        // Upper diagonal
        if (j < m - 2) {
            Q_i[j][j + 1] = (1.0 - theta) * current_c / (delta_x * delta_x);
        }
    }
    return Matrix(Q_i);
}

Matrix PDEPricer::getVi(int timestep_i, double delta_x){
    // Getting the params
    std::vector<double> b_i = B.data[timestep_i];
    std::vector<double> c_i = C.data[timestep_i];
    std::vector<double> d_i = D.data[timestep_i];

    // Initialize V_i as a (m-1) rows x 1 col vector
    std::vector<std::vector<double>> V_i(m - 1, std::vector<double>(1, 0.0));

    // Getting the boundary conditions
    double u_i_xm = bound_cond[timestep_i][0]; double u_i1_xm = bound_cond[timestep_i + 1][0];
    double u_i_xM = bound_cond[timestep_i][1]; double u_i1_xM = bound_cond[timestep_i + 1][1];

    // Fill the vector
    for (int j = 0; j < m - 1; j++) {
        double current_c = c_i[j + 1]; double current_b = b_i[j + 1];
        V_i[j][0] = d_i[j + 1]; // Main element
        // Lower bound condition
        if (j==0){
            V_i[j][0] += (-current_b / (2.0 * delta_x) + theta * current_c / (delta_x * delta_x)) * u_i_xm +
                    (1.0 - theta) * current_c * u_i1_xm / (delta_x * delta_x);
        }
        // Upper bound condition
        if (j==m-2){
            V_i[j][0] += (current_b / (2.0 * delta_x) + theta * current_c / (delta_x * delta_x)) * u_i_xM +
                    (1.0 - theta) * current_c * u_i1_xM / (delta_x * delta_x);
        }
    }
    return Matrix(V_i);
}

std::vector<Matrix> PDEPricer::getPricePaths(double delta_x) {
    // Initializing the vector of vectors containing all the U_i
    std::vector<Matrix> U(n, Matrix());

    // Initializing the vector U using terminal condition - removing the extreme points to get a matrix of size (m-1)²
    std::vector<double> term_cond_interior = term_cond;
    term_cond_interior.erase(term_cond_interior.begin());
    term_cond_interior.erase(term_cond_interior.end() - 1);

    std::vector<std::vector<double>> mat_term_cond_interior(m - 1, std::vector<double>(1, 0.0));
    for (int i = 0; i < m - 1; ++i) {
        mat_term_cond_interior[i][0] = term_cond_interior[i];
    }
    Matrix U_i1 = Matrix(mat_term_cond_interior);

    // Instantiating a matrix of 0 to perform substraction
    Matrix Null_matrix(std::vector<std::vector<double>>(m - 1, std::vector<double>(m - 1, 0.0)));

    for (int i = (n - 1); i >= 0; i--) { // Going through time grid backward
        Matrix Pi = getPi(i, delta_x);
        Matrix Q_i = getQi(i, delta_x);
        Matrix V_i = getVi(i, delta_x);
        Matrix P_inv = Pi.inverseLU();
        Matrix P_inv_neg = Null_matrix - P_inv;
        Matrix interm = (Q_i * U_i1) + V_i;
        Matrix U_i = P_inv_neg * interm;
        U_i1.setMatrixData(U_i);
        U[i] = U_i;
    }
    return U;
}




#include <iostream>

double PDEPricer::getOptionPrice(double initialSpot) {
    // Check if space_grid is not empty and has at least two elements
    if (space_grid.size() < 2) {
        throw std::runtime_error("space_grid is too small.");
    }

    double delta_x = space_grid[1] - space_grid[0];

    std::vector<Matrix> optionPrices = getPricePaths(delta_x);

    if (optionPrices.size() != static_cast<std::size_t>(this->n)) {
        throw std::runtime_error("Error in the resolution path of the PDE.");
    } // Check if we have all the paths from t = 0 to n-1; t=T is known by terminal condition

    double min_elem = *std::min_element(space_grid.begin(), space_grid.end());
    double max_elem = *std::max_element(space_grid.begin(), space_grid.end());

    if ((initialSpot < min_elem) || (initialSpot > max_elem)) {
        throw std::runtime_error("Please adjust the space grid, boundaries are not valid.");
    }

    std::vector<std::vector<double>> vecPrice = optionPrices[0].data;

    // Finding the index in the space grid of the elem closest to S0
    int closest_idx = 0;
    for (std::size_t i = 0; i < space_grid.size(); i++) {
        double difference = space_grid[i] - initialSpot;
        if (difference > 0) { /* take the first spot price higher than the spot */
            closest_idx = i;
            break;
        }
    }

    // if possible, we do a linear interpolation between the two option price, otherwise just return the closest price
    double fd_price = vecPrice[closest_idx][0];
    double threshold = 1e-6;
    if ((closest_idx > 0) && (std::abs(space_grid[closest_idx] - initialSpot) > threshold)) {
        double d = (initialSpot - space_grid[closest_idx - 1]) / delta_x;
        fd_price = (1 - d) * vecPrice[closest_idx - 1][0] + d * vecPrice[closest_idx][0];
    }

    return fd_price;
}
