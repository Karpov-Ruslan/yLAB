#include <fstream>
#include <iostream>
#include "geometry.hpp"

namespace krvlib {

void test(int exp, const triangle& triangle_1, const triangle& triangle_2) {
    std::cout << "Expected " << exp << ": " << triangle::intersect(triangle_1, triangle_2) << std::endl;
}

}

using namespace krvlib;

int main() {
    
    std::ifstream fout("../tests.txt");

    if (!fout.is_open()) {
        std::cout << "file \"../tests.txt\" did not open" << std::endl;
        return 1;
    }

    int num_of_tests;
    fout >> num_of_tests;

    for (int i = 0; i < num_of_tests; i++) {
        int expected;
        double point_1[9];
        double point_2[9];
        fout >> expected;
        for (int j = 0; j < 9; j++) {
            fout >> point_1[j];
        }
        for (int j = 0; j < 9; j++) {
            fout >> point_2[j];
        }
        test(expected, triangle(vector(point_1[0], point_1[1], point_1[2]), vector(point_1[3], point_1[4], point_1[5]), vector(point_1[6], point_1[7], point_1[8])),
                       triangle(vector(point_2[0], point_2[1], point_2[2]), vector(point_2[3], point_2[4], point_2[5]), vector(point_2[6], point_2[7], point_2[8])));
        
    }

    fout.close();
    return 0;
}