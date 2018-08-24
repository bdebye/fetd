
#ifndef EIGEN_H_
#define EIGEN_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <iostream>

using namespace Eigen;

inline double norm(const Vector3d &vec) {
    return vec.norm();
    //return sqrt(pow(vec(0), 2) + pow(vec(1), 2) + pow(vec(2), 2));
}

inline double dot(const Vector3d &a, const Vector3d &b) {
    return a.dot(b);
}

inline Vector3d cross(const Vector3d &a, const Vector3d &b) {
    return a.cross(b);
}

template <typename T> inline void print(const T &value) {
    std::cout << value << std::endl;
}

#endif /* EIGEN_H_ */
