double MCEuCallPricer(double S0, double K, double T, double r, double q, double sigma, int numSimulations);

class QMCEuCallPricer {
public:
    // Initializing constructor
    QMCEuCallPricer(double S0, double K, double T, double r, double q, double sigma, int numSimulations);

    // Member Functions
    void shuffleVector(std::vector<double>& vec);
    void initialize_sobol_sequence();
    void initialize_quasi_random_normal_variables();
    double get_qmc_price();

private:
    double S0; double K; double T; double r; double q; double sigma;
    int numSimulations;
    std::vector<double> sobolSequence; std::vector<double> normalSequence;
};