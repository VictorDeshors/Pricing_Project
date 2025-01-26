# Pricing of a vanilla European option

C++ Project I did after _MAP552 Stochastic Calculus_ on the pricing of a vanilla European option.

The case of a simple European call option enables us to introduce the following concepts: 
1. Black-Scholes (B.S.) formula
2. Finite difference method of the B.S. partial differential equation (P.D.E.)
3. Monte Carlo & Quasi-Monte Carlo methods using Sobol suites

Each method is explained in detail in [Theory.pdf](Theory.pdf).
The plots are done in python in [plot.ipynb](plot.ipynb).

## Structure of the C++ code
- [matrix.cpp](matrix.cpp): implements a Matrix class performing basic operations and computing its inverse using LU decomposition.
- [pde_pricer.cpp](pde_pricer.cpp): implements a pricer solving in a backward way the P.D.E.
- [simulations.cpp](simulations.cpp): implements a Monte Carlo pricer as well as a Quasi-Monte carlo pricing class.

## Building and Running the C++ Code

To build and run the C++ code, follow these steps:

1. Create a build directory and navigate into it:
   ```sh
   mkdir build
   cd build
   ```

2. Run CMake to configure the project:
    ```sh
    cmake ..
    ```

3. Build the project:
    ```sh
    make
    ```

4. Run the executable:
    ```sh
    ./PDE_Victor_DESHORS
    ```
