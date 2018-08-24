

#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_

#include <vector>
#include <array>

#include "element.h"
#include "boundary.h"

extern std::vector<element> tetrah_element;
extern std::vector<std::array<int, 2>> edge_basis;
extern std::vector<bdy_element> boundary_element;
extern VectorXd basis_coef;

extern SparseMatrix<double> M;
extern SparseMatrix<double> P;
extern SparseMatrix<double> W;
extern SparseMatrix<double> K;

void parse_edge_basis();

void generate_elements();

void generate_boundary_elements();

void print_edge_basis();

int find_edge(int node_a, int node_b);

void assembly();

Vector3d t(int I);

void solve_eigenmode();

void newmark(int n);

#endif /* ASSEMBLY_H_ */
