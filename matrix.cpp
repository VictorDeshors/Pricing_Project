#include "matrix.h"
#include <iostream>
#include <vector>
#include <cmath>
//#include <Eigen/Dense>

// Default constructor - empty 2*2 matrix
Matrix::Matrix(): data(2, std::vector<double>(2, 0.0)) {}

// Parameterized constructor
Matrix::Matrix(const std::vector<std::vector<double>>& vect): data(vect) {}

void Matrix::setMatrixData(const Matrix& other_mat) {
    this->data = other_mat.data;
}
void Matrix::addRow(const std::vector<double>& newRow, int newRowIndex) {
    if (newRowIndex >= 0 && static_cast<std::size_t>(newRowIndex) <= data.size()) { // Ensure n is within valid bounds
        data.insert(data.begin() + newRowIndex, newRow);
    } else {
        std::cerr << "Index out of bounds!" << std::endl;
    }
}

double Matrix::getMatrixDeterminant(){
    return computeDeterminant(data);
}

double Matrix::computeDeterminant(std::vector<std::vector<double>> vect){
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } // Check if the matrix is quadratic
    int dimension = vect.size();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect[0][0];
    }

    // Formula for 2x2-matrix
    if(dimension == 2) {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }
    double result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {
        // Submatrix
        std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double> (dimension - 1));
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect[m-1][z] = data[m][n];
                    z++;
                }
            }
        }
        //recursive call
        result = result + sign * data[0][i] * computeDeterminant(subVect);
        sign = -sign;
    }

    return result;
}

std::vector<std::vector<double>> Matrix::getTranspose(const std::vector<std::vector<double>> matrix1) {
    // Creating a matrix of the same size
    std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double> (matrix1.size()));

    // Filling solution-matrix
    for(size_t i = 0; i < matrix1.size(); i++) {
        for(size_t j = 0; j < matrix1[0].size(); j++) {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

std::vector<std::vector<double>> Matrix::getCofactor(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    }

    std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
    std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * computeDeterminant(subVect);
        }
    }
    return solution;
}

Matrix Matrix::getInverseMatrix() {
    double det = computeDeterminant(data);
    if(det == 0) {
        throw std::runtime_error("Determinant is 0");
    } // if matrix is not invertible, throw error
    double d = 1.0/computeDeterminant(data);
    std::vector<std::vector<double>> solution(data.size(), std::vector<double> (data.size()));

    for(size_t i = 0; i < data.size(); i++) {
        for(size_t j = 0; j < data.size(); j++) {
            solution[i][j] = data[i][j];
        }
    }
    solution = getTranspose(getCofactor(solution));
    // Scaling by the inverse determinant factor
    for(size_t i = 0; i < data.size(); i++) {
        for(size_t j = 0; j < data.size(); j++) {
            solution[i][j] *= d;
        }
    }
    return Matrix(solution);
}

bool Matrix::luDecomposition(const std::vector<std::vector<double>>& matrix,
                     std::vector<std::vector<double>>& L,
                     std::vector<std::vector<double>>& U) {
    size_t n = matrix.size();

    // Initialize L and U
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        // Upper triangular matrix U
        for (size_t k = i; k < n; ++k) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = matrix[i][k] - sum;
        }

        // Lower triangular matrix L
        for (size_t k = i; k < n; ++k) {
            if (i == k) {
                L[i][i] = 1.0; // Diagonal of L is 1
            } else {
                double sum = 0.0;
                for (size_t j = 0; j < i; ++j) {
                    sum += L[k][j] * U[j][i];
                }
                if (std::fabs(U[i][i]) < 1e-9) {
                    return false; // Singular matrix
                }
                L[k][i] = (matrix[k][i] - sum) / U[i][i];
            }
        }
    }

    return true;
}

double Matrix::determinantLU() {
    size_t n = data.size();
    std::vector<std::vector<double>> L, U;

    if (!luDecomposition(data, L, U)) {
        return 0.0; // Singular matrix
    }

    double det = 1.0;
    for (size_t i = 0; i < n; ++i) {
        det *= U[i][i];
    }

    return det;
}

Matrix Matrix::inverseLU() {
    size_t n = data.size();
    std::vector<std::vector<double>> L, U;

    if (!luDecomposition(data, L, U)) {
        throw std::runtime_error("Matrix is singular, inverse does not exist.");
    }
    // Initialize identity matrix
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        I[i][i] = 1.0;
    }
    // Forward substitution to solve L * Y = I
    std::vector<std::vector<double>> Y(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L[i][k] * Y[k][j];
            }
            Y[i][j] = I[i][j] - sum;
        }
    }
    // Back substitution to solve U * X = Y
    std::vector<std::vector<double>> X(n, std::vector<double>(n, 0.0));
    for (int i = n - 1; i >= 0; --i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = i + 1; k < n; ++k) {
                sum += U[i][k] * X[k][j];
            }
            X[i][j] = (Y[i][j] - sum) / U[i][i];
        }
    }
    Matrix InvMat = Matrix(X);
    return InvMat;
}

//// Convert std::vector<std::vector<double>> to Eigen::MatrixXd
//Eigen::MatrixXd Matrix::vectorToEigen(const std::vector<std::vector<double>>& vec) {
//    size_t rows = vec.size();
//    size_t cols = vec[0].size();
//    Eigen::MatrixXd mat(rows, cols);
//    for (size_t i = 0; i < rows; ++i) {
//        for (size_t j = 0; j < cols; ++j) {
//            mat(i, j) = vec[i][j];
//        }
//    }
//    return mat;
//}
//
//// Convert Eigen::MatrixXd to std::vector<std::vector<double>>
//std::vector<std::vector<double>> Matrix::eigenToVector(const Eigen::MatrixXd& mat) {
//    size_t rows = mat.rows();
//    size_t cols = mat.cols();
//    std::vector<std::vector<double>> vec(rows, std::vector<double>(cols));
//    for (size_t i = 0; i < rows; ++i) {
//        for (size_t j = 0; j < cols; ++j) {
//            vec[i][j] = mat(i, j);
//        }
//    }
//    return vec;
//}
//
//// Function to compute the inverse of a matrix
//Matrix Matrix::computeInverseEigen() {
//    Eigen::MatrixXd mat = vectorToEigen(data);
//
//    if (mat.rows() != mat.cols()) {
//        throw std::invalid_argument("Matrix must be square to compute the inverse.");
//    }
//
//    Eigen::MatrixXd invMat = mat.inverse();
//    return Matrix(eigenToVector(invMat));
//}

Matrix Matrix::operator-(const Matrix& other) const {
    // Check nb of rows and cols
    if (data.size() != other.data.size()) {
        throw std::runtime_error("Matrices don't have the same number of rows.");
    }
    if (data[0].size() != other.data[0].size()) {
        throw std::runtime_error("Matrices don't have the same number of columns.");
    }
    std::vector<std::vector<double>> diff_vect(data.size(), std::vector<double>(data[0].size()));
    for (int i = 0; i < static_cast<int>(data.size()); i++) {
        for (int j = 0; j < static_cast<int>(data[0].size()); j++) {
            diff_vect[i][j] = data[i][j] - other.data[i][j];
        }
    }
    return Matrix(diff_vect);
}

Matrix Matrix::operator+(const Matrix& other) const {
    // Check nb of rows and cols
    if (data.size() != other.data.size()) {
        throw std::runtime_error("Matrices don't have the same number of rows.");
    }
    if (data[0].size() != other.data[0].size()) {
        throw std::runtime_error("Matrices don't have the same number of columns.");
    }
    std::vector<std::vector<double>> sum_vect(data.size(), std::vector<double>(data[0].size()));
    for (int i = 0; i < static_cast<int>(data.size()); i++) {
        for (int j = 0; j < static_cast<int>(data[0].size()); j++) {
            sum_vect[i][j] = data[i][j] + other.data[i][j];
        }
    }
    return Matrix(sum_vect);
}

Matrix Matrix::operator*(const Matrix& other) const {
    // Ensure matrices are non-empty
    if (data.empty() || other.data.empty() || other.data[0].empty()) {
        throw std::runtime_error("Matrices cannot be multiplied: one or both matrices are empty.");
    }

    // Check dimension compatibility
    if (data[0].size() != other.data.size()) {
        throw std::runtime_error("Matrices cannot be multiplied: number of columns in the first matrix != number of rows in the second matrix.");
    }

    // Initialize the result matrix
    size_t rows = data.size();
    size_t cols = other.data[0].size();
    size_t common_dim = other.data.size();
    std::vector<std::vector<double>> prod_vect(rows, std::vector<double>(cols, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            for (size_t k = 0; k < common_dim; ++k) {
                prod_vect[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    return Matrix(prod_vect);
}