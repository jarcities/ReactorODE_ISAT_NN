#include <iostream>
#include <array>
#include <iomanip>  // for std::setprecision
#include "fnn_codejenn.h"  // Your generated header with cnn2(...) definition

using Scalar = double;

int main() {

    std::array<Scalar, 11> input = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    

    // Pass the input to your generated CNN function
    auto output = myfnn<Scalar>(input);

    // Print the results with high precision
    std::cout << std::scientific << std::setprecision(15);  // Set precision and scientific notation
    std::cout << "Output:\n";  // Print each value on a new line
    for(const auto& val : output) {
        std::cout << val << '\n';
    }
    std::cout << std::endl;

    return 0;
}

/*
Compile and run:
clang++ -std=c++23 -Wall -O3 -march=native -o test test.cpp
./test
*/