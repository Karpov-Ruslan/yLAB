#include <iostream>
#include <vector>
#include "geometry.hpp"

using namespace krvlib;

int main() {

    std::vector<triangle> arr;
    int num;
    std::cin >> num;
    bool intersects[num] = {};

    for (int i = 0; i < num; i++) {
        double point[9];
        for (int j = 0; j < 9; j++) {
            std::cin >> point[j];
        }
        arr.push_back(triangle(vector(point[0], point[1], point[2]), vector(point[3], point[4], point[5]), vector(point[6], point[7], point[8])));
    }


    for (int i = 0; i < num; i++) {
        if (intersects[i]) {
            continue;
        }
        for (int j = 0; j < num; j++) {
            if (i == j) {continue;}
            if (triangle::intersect(arr[i], arr[j])) {
                intersects[i] = intersects[j] = true;
            }
        }
    }
    
    int count = 0;
    for (int i = 0; i < num; i++) {
        if (intersects[i]) {
            std::cout << i + 1 << " ";
        }
    }

    return 0;
}