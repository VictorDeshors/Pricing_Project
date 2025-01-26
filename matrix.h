#include <vector>
//#include <Eigen/Dense>

class Matrix {
public:
    std::vector<std::vector<double>> data;
    // Constructors
    Matrix(); // Default constructor
    explicit Matrix(const std::vector<std::vector<double>>& vect); // Initializing constructor

    // Member Functions
    void setMatrixData(const Matrix& other_mat);
    void addRow(const std::vector<double>& newRow, int newRowIndex);
    double computeDeterminant(std::vector<std::vector<double>> vect);
    double getMatrixDeterminant();
    std::vector<std::vector<double>> getTranspose(std::vector<std::vector<double>> vect);
    std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>> vect);
    Matrix getInverseMatrix();

    bool luDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U);
    double determinantLU();
    Matrix inverseLU();

//    Eigen::MatrixXd vectorToEigen(const std::vector<std::vector<double>>& vec);
//    std::vector<std::vector<double>> eigenToVector(const Eigen::MatrixXd& mat);
    Matrix computeInverseEigen();


        // Overloaded Operators
    Matrix operator-(const Matrix &other) const; // Difference of two matrix objects
    Matrix operator+(const Matrix &other) const; // Sums two matrix objects
    Matrix operator*(const Matrix &other) const; // Multiply two matrix objects

private:
};