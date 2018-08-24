
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "eigen.h"
#include <array>
#include <vector>

extern const int BDY_GROUP_TAG;

struct bdy_element {
    
    int bdy_element_index;
    
    std::array<int, 3> basis_global_index;
    
    double area;
    
    double epsilon;
    double sigma;
    double mu;
    
    std::array<Vector2d, 3> vertex;
    Matrix3d coef_matrix;
    
    static bdy_element creat_boundary_element(int n);
    
    double a(int i) {
        return coef_matrix(i-1, 2);
    }
    
    double b(int i) {
        return coef_matrix(i-1, 0);
    }
    
    double c(int i) {
        return coef_matrix(i-1, 1);
    }
    
    double p(int i, int j) {
        return b(i) * b(j) + c(i) * c(j);
    }
    
    double b_1_Jin() {
        return (vertex[1](1) - vertex[2](1)) / (2 * area);
    }
    
    double l(int e);
    double F(int i, int j);
    double W(int i, int j);
    
    double Y() {
        return sqrt(epsilon / mu);
    }
    
    int simplex_index(int e, int i);
    
};

//extern std::vector<bdy_element> bdy_element;

#endif /* BOUNDARY_H_ */
