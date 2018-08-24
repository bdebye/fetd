

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "eigen.h"
#include "gmsh.h"

#include <array>

struct element {
    
    Matrix4d coef_matrix;
    
    double volumn;

    
    std::array<int, 6> basis_global_index;
    int element_index;
    
    
    double epsilon;
    double sigma;
    double mu;
    
    double a(int i) { return coef_matrix(i-1, 3); }
    double b(int i) { return coef_matrix(i-1, 0); }
    double c(int i) { return coef_matrix(i-1, 1); }
    double d(int i) { return coef_matrix(i-1, 2); }
    
    Vector3d u(int i, int j) {
        Vector3d upsi;
        upsi(0) = c(i) * d(j) - c(j) * d(i);
        upsi(1) = b(j) * d(i) - b(i) * d(j);
        upsi(2) = b(i) * c(j) - b(j) * c(i);
        return upsi;
    }
    
    double p(int i, int j) {
        return b(i) * b(j) + c(i) * c(j) + d(i) * d(j);
    }
    
    Vector3d t(int e);
    double l(int e);
    
    int simplex_index(int e, int i);
    
    double E(int i, int j);
    double F(int i, int j);
    
    double K(int i, int j);
    double M(int i, int j);
    double P(int i, int j);
    
    static element create_element(int tet_n);
    
    Vector3d center_coord();
    Vector3d center_field();
    
};

#endif /* ELEMENT_H_ */
