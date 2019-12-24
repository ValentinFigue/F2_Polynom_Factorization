cmake_minimum_required(VERSION 3.1)  # CMake version check
project(F2_Polynom_Factorization)               

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(F2_Polynom_Factorization ${SOURCE_FILES})
