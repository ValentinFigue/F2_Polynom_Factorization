# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project(F2_Polynom_Factorization)               

# include
include_directories( ./include/F2_Polynom_Factorization ./src)

# Add main.cpp file of project root directory as source file
add_executable(F2_Polynom_Factorization main.cpp ./src/Arithmetic.cpp ./src/F2_Polynom.cpp ./src/utils.cpp ./include/F2_Polynom_Factorization/utils.hpp ./include/F2_Polynom_Factorization/F2_Polynom.hpp ./include/F2_Polynom_Factorization/Arithmetic.hpp)

