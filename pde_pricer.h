#include <vector>
#include "matrix.h"

class PDEPricer {
public:
    PDEPricer(
            double T, int n, int m, std::vector<double> space_grid, const std::vector<std::vector<double>>& bound_cond,
            const std::vector<double>& term_cond, const Matrix& A, const Matrix& B, const Matrix& C, const Matrix& D,
            double theta
            ); // Initializing constructor

    // Member Functions
    Matrix getPi(int timestep_i, double delta_x);
    Matrix getQi(int timestep_i, double delta_x);
    Matrix getVi(int timestep_i, double delta_x);
    std::vector<Matrix> getPricePaths(double delta_x);
    double getOptionPrice(double initialSpot);


private:
    double T; // Expiry of the option; used to create the time grid
    int n; // nb of steps to discretize the time grid
    int m; // nb of steps to discretize the space grid
    std::vector<double> space_grid;
    std::vector<std::vector<double>> bound_cond; // should be of size (n, 2)
    std::vector<double> term_cond;  // should be of size (n, 1)
    Matrix A, B, C, D;
    double theta; // smoothing param
};