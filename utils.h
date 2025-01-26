double blackScholesPrice(double S, double K, double T, double r, double sigma, std::string optionType);
std::vector<double> getGrid(double init, double step, int length);
void printVector(std::vector<double> vec);