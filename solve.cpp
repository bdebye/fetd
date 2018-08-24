

#include "solve.h"
#include "eigen.h"
#include "gmsh.h"
#include "element.h"
#include "dataio.h"
#include "feed.h"
#include "waveform.h"

#include <iostream>

std::vector<std::array<int, 2>> edge_basis;
std::vector<element> tetrah_element;
std::vector<bdy_element> boundary_element;
VectorXd basis_coef;

int find_edge(int node_a, int node_b) {
    for(int i = 0; i < edge_basis.size(); i++) {
        if((edge_basis[i][0] == node_a && edge_basis[i][1] == node_b)
           || (edge_basis[i][1] == node_a && edge_basis[i][0] == node_b)) {
            return i;
        }
    }
    return -1;
}

bool edge_already_basis(int node_a, int node_b) {
    return find_edge(node_a, node_b) >= 0;
}

void push_edge_basis(int node_a, int node_b) {
    if(!edge_already_basis(node_a, node_b)) {
        std::array<int, 2> basis;
        basis[0] = node_a;
        basis[1] = node_b;
        edge_basis.push_back(basis);
    }
}

void parse_edge_basis() {
    for(int n = 0; n < mesh.tetrah_num; n++) {
        const msh_tetrah &tetrah = mesh.tetrah[n];
        push_edge_basis(tetrah.node_number_list[0], tetrah.node_number_list[1]);
        push_edge_basis(tetrah.node_number_list[0], tetrah.node_number_list[2]);
        push_edge_basis(tetrah.node_number_list[0], tetrah.node_number_list[3]);
        push_edge_basis(tetrah.node_number_list[1], tetrah.node_number_list[2]);
        push_edge_basis(tetrah.node_number_list[1], tetrah.node_number_list[3]);
        push_edge_basis(tetrah.node_number_list[2], tetrah.node_number_list[3]);
    }
    std::cout << "Parsed the edges in the mesh, " << edge_basis.size() << " edges found." << std::endl;
}

void print_edge_basis() {
    for(int i = 0; i < edge_basis.size(); i++) {
        std::cout << i << ": " << edge_basis[i][0] << "<----->" << edge_basis[i][1] << std::endl;
    }
}

double edge_length(int n) {
    int node_a = edge_basis[n][0];
    int node_b = edge_basis[n][1];
    return norm(get_node(node_a) - get_node(node_b));
}


void generate_elements() {
    if(edge_basis.size() == 0) {
        std::cout << "The edges haven't been parsed..." << std::endl;
    }
    tetrah_element.clear();
    tetrah_element.reserve(mesh.tetrah_num);
    for(int i = 0; i < mesh.tetrah_num; i++) {
        tetrah_element.push_back(element::create_element(i));
    }
    generate_boundary_elements();
}

void generate_boundary_elements() {
    if(edge_basis.size() == 0) {
        std::cout << "The edges haven't been parsed..." << std::endl;
    }
    boundary_element.clear();
    boundary_element.reserve(mesh.triang_num);
    for(int i = 0; i < mesh.triang_num; i++) {
        boundary_element.push_back(bdy_element::creat_boundary_element(i));
    }
}

Vector3d t(int I) {
    int node_a = edge_basis[I][0];
    int node_b = edge_basis[I][1];
    Vector3d edge_v = get_node(node_a) - get_node(node_b);
    return edge_v / norm(edge_v);
}


SparseMatrix<double> M;
SparseMatrix<double> P;
SparseMatrix<double> W;
SparseMatrix<double> K;
VectorXd h;

/*
void assembly() {
    int N = (int) edge_basis.size();
    M.resize(N, N);
    P.resize(N, N);
    W.resize(N, N);
    K.resize(N, N);
    h.resize(N);
    std::cout << "Assembly the equations..." << std::endl;
    for(int e = 0; e < tetrah_element.size(); e++) {
        for(int i = 1; i <= 6; i++) {
            for(int j = 1; j <= 6; j++) {
                element &elem = tetrah_element[e];
                int I = elem.basis_global_index[i-1];
                int J = elem.basis_global_index[j-1];
                M.coeffRef(I, J) += elem.M(i, j);
                K.coeffRef(I, J) += elem.K(i, j);
                P.coeffRef(I, J) += elem.P(i, j);
            }
        }
    }
    //std::cout << M << std::endl;
}
 */

void assembly_boundary() {
    int N = (int) edge_basis.size();
    W.resize(N, N);
    std::vector<Eigen::Triplet<double>> W_triple_list;
    for(int e = 0; e < boundary_element.size(); e++) {
        for(int i = 1; i <= 3; i++) {
            for(int j = 1; j <= 3; j++) {
                bdy_element &bdy_elem = boundary_element[e];
                int I = bdy_elem.basis_global_index[i-1];
                int J = bdy_elem.basis_global_index[j-1];
                W_triple_list.push_back(Eigen::Triplet<double>(I, J, bdy_elem.W(i, j)));
            }
        }
    }
    W.setFromTriplets(W_triple_list.begin(), W_triple_list.end());
}

void assembly() {
    int N = (int) edge_basis.size();
    M.resize(N, N);
    P.resize(N, N);
    W.resize(N, N);
    K.resize(N, N);
    h.resize(N);
    
    std::cout << "Assembly the equations..." << std::endl;
    std::vector<Eigen::Triplet<double>> M_triple_list;
    std::vector<Eigen::Triplet<double>> K_triple_list;
    std::vector<Eigen::Triplet<double>> P_triple_list;
    
    for(int e = 0; e < tetrah_element.size(); e++) {
        for(int i = 1; i <= 6; i++) {
            for(int j = 1; j <= 6; j++) {
                element &elem = tetrah_element[e];
                int I = elem.basis_global_index[i-1];
                int J = elem.basis_global_index[j-1];
                M_triple_list.push_back(Eigen::Triplet<double>(I, J, elem.M(i, j)));
                K_triple_list.push_back(Eigen::Triplet<double>(I, J, elem.K(i, j)));
                P_triple_list.push_back(Eigen::Triplet<double>(I, J, elem.P(i, j)));
            }
        }
    }
    M.setFromTriplets(M_triple_list.begin(), M_triple_list.end());
    K.setFromTriplets(K_triple_list.begin(), K_triple_list.end());
    P.setFromTriplets(P_triple_list.begin(), P_triple_list.end());
    
    assembly_boundary();
    //std::cout << M << std::endl;
}

void solve_eigenmode() {
    GeneralizedEigenSolver<MatrixXd> eigen_solver;
    std::cout << "Compute the eigenvalue system, the dimension of the problem is "
                << edge_basis.size() << "." << std::endl;
    eigen_solver.compute(K, M);
    std::cout << eigen_solver.eigenvalues() << std::endl;
}


const double dt = 1e-11;

const double gamma_ = 0.5;
const double beta_ = 0.25;
const double pi_ = M_PI;

const double f = 1e8;

#define DEFAULT_TAU 2e-9

double excitation(double t) {
    //return gauss_dd(t - 5 * DEFAULT_TAU, DEFAULT_TAU);
    return cos(2 * pi_ * f * t);
}

MatrixXd sparse_P_inverse(const SparseMatrix<double> A) {
    int N = (int) edge_basis.size();
    SimplicialLLT<SparseMatrix<double>> solver;
    solver.compute(A);
    MatrixXd I(N, N);
    I.setIdentity();
    return solver.solve(I);
}

void newmark(int n) {
    int N_edge = (int) edge_basis.size();
    SparseMatrix<double> C = P + W;
    
    SparseMatrix<double> P_ = M + gamma_ * dt * C + beta_ * pow(dt, 2) * K;
    SparseMatrix<double> Q_ = -2.0 * M + (1.0 - 2.0 * gamma_) * dt * C +
                            (0.5 + gamma_ - 2 * beta_) * pow(dt, 2) * K;
    SparseMatrix<double> R_ = M + (gamma_ - 1) * dt * C + (0.5 - gamma_ + beta_) * pow(dt, 2) * K;
    std::cout << "Solving the inverse of P matrix..." << std::endl;
    MatrixXd INVP_ = sparse_P_inverse(P_);
    
    VectorXd h1(N_edge);
    VectorXd h2(N_edge);
    VectorXd h3(N_edge);
    
    VectorXd ph1(N_edge);
    VectorXd ph2(N_edge);
    VectorXd ph3(N_edge);
    VectorXd ph_temp(N_edge);
    
    for(int i = 0; i < N_edge; i++) {
        h1(i) = 0;
        h2(i) = 0;
        h3(i) = 0;
        ph1(i) = 0;
        ph2(i) = 0;
        ph3(i) = 0;
        ph_temp(i) = 0;
    }
    
    struct signal receive;
    receive.dt = dt;
    const int SOURCE = find_feed_index();
    //std::cout << SOURCE << std::endl;
    double t = 0;
    for(int i = 0; i <= n; i++) {
        t = t + dt;
        std::cout << i << ": " << excitation(t) << std::endl;
        
        h3(SOURCE) = - edge_length(SOURCE) * excitation(t);
        ph_temp = - Q_ * ph2 - R_ *  ph1
        + pow(dt, 2) * (beta_ * h3 + (0.5 + gamma_ - 2 * beta_)
                        * h2 + (0.5 - gamma_ + beta_) * h1);
        ph3 = INVP_ * ph_temp;
        h1 = h2;
        h2 = h3;
        ph1 = ph2;
        ph2 = ph3;
        basis_coef = ph3;
        
        receive.push_back(tetrah_element[1637].center_field()(2));
        //receive.push_back(tetrah_element[2891].center_field()(2));
        
        if(i % 500 == 0) {
            export_scalar_field("scalar_field_" + integer_to_string(i) + ".pos");
            export_vector_field("vector_field_" + integer_to_string(i) + ".pos");
        }
    }
    receive.save("signal.txt");
}

