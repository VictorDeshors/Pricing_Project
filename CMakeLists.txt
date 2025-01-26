cmake_minimum_required(VERSION 3.30)
project(PDE_Victor_DESHORS)
set(CMAKE_CXX_STANDARD 20)

add_executable(PDE_Victor_DESHORS main.cpp
        matrix.cpp
        matrix.h
        utils.cpp
        utils.h
        pde_pricer.cpp
        pde_pricer.h
        simulations.cpp
        simulations.h)
