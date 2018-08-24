
#include "element.h"
#include "gmsh.h"
#include "solve.h"


int element::simplex_index(int e, int i) {
    const static int index[6][2] = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};
    return index[e-1][i-1];
}

element element::create_element(int tet_n) {
    element elem;
    msh_tetrah tetrah = mesh.tetrah[tet_n];
    std::array<Vector3d, 4> vertex;
    elem.element_index = tet_n;
    for(int i = 0; i < 4; i++) {
        vertex[i] = mesh.node[tetrah.node_number_list[i]];
    }
    if(medium_conf.find(tetrah.physical_tag) != medium_conf.end()) {
        elem.epsilon = medium_conf[tetrah.physical_tag].epsilon;
        elem.sigma = medium_conf[tetrah.physical_tag].sigma;
        elem.mu = medium_conf[tetrah.physical_tag].mu;
    }
    else {
        elem.epsilon = 8.854e-012;
        elem.sigma = 0;
        elem.mu = 1.257e-006;
    }
    
    elem.basis_global_index[0] = find_edge(tetrah.node_number_list[0], tetrah.node_number_list[1]);
    elem.basis_global_index[1] = find_edge(tetrah.node_number_list[0], tetrah.node_number_list[2]);
    elem.basis_global_index[2] = find_edge(tetrah.node_number_list[0], tetrah.node_number_list[3]);
    elem.basis_global_index[3] = find_edge(tetrah.node_number_list[1], tetrah.node_number_list[2]);
    elem.basis_global_index[4] = find_edge(tetrah.node_number_list[1], tetrah.node_number_list[3]);
    elem.basis_global_index[5] = find_edge(tetrah.node_number_list[2], tetrah.node_number_list[3]);
    
    Matrix4d coord_matrix;
    for(int j = 0; j < 4; j++) {
        coord_matrix(0, j) = vertex[j](0);
        coord_matrix(1, j) = vertex[j](1);
        coord_matrix(2, j) = vertex[j](2);
        coord_matrix(3, j) = 1;
    }
    // std::cout << coord_matrix << std::endl;
    elem.coef_matrix = coord_matrix.inverse();
    elem.volumn = fabs(coord_matrix.determinant()) / 6.0;
    
    return elem;
}

double element::E(int i, int j) {
    int i1 = simplex_index(i, 1);
    int i2 = simplex_index(i, 2);
    int j1 = simplex_index(j, 1);
    int j2 = simplex_index(j, 2);
    return 4.0 * volumn * l(i) * l(j) * dot(u(i1, i2), u(j1, j2));
}

double element::F(int i, int j) {
    const static Matrix4d M = (Matrix4d::Ones() + Matrix4d::Identity()) / 20.0;
    int i1 = simplex_index(i, 1);
    int i2 = simplex_index(i, 2);
    int j1 = simplex_index(j, 1);
    int j2 = simplex_index(j, 2);
    return volumn * l(i) * l(j)
        * (p(i2, j2) * M(i1-1, j1-1)
         - p(i2, j1) * M(i1-1, j2-1)
         - p(i1, j2) * M(i2-1, j1-1)
         + p(i1, j1) * M(i2-1, j2-1));
}

double element::l(int i) {
    int node_a = edge_basis[this->basis_global_index[i-1]][0];
    int node_b = edge_basis[this->basis_global_index[i-1]][1];
    return norm(get_node(node_a) - get_node(node_b));
}

Vector3d element::t(int e) {
    int node_a = edge_basis[basis_global_index[e-1]][0];
    int node_b = edge_basis[basis_global_index[e-1]][1];
    Vector3d basis_v = get_node(node_a) - get_node(node_b);
    return basis_v / norm(basis_v);
}

double element::K(int i, int j) {
    return this->E(i, j) / this->mu;
}

double element::M(int i, int j) {
    return this->F(i, j) * this->epsilon;
}

double element::P(int i, int j) {
    return this->F(i, j) * this->sigma;
}

Vector3d element::center_coord() {
    Vector3d center = {0, 0, 0};
    for(int i = 0; i < 4; i++) {
        center += get_node(mesh.tetrah[this->element_index].node_number_list[i]);
    }
    return center / 4.0;
}

Vector3d element::center_field() {
    Vector3d field = {0, 0, 0};
    Vector3d center = this->center_coord();
    double x = center(0);
    double y = center(1);
    double z = center(2);
    for(int i = 1; i <= 6; i++) {
        int i1 = this->simplex_index(i, 1);
        int i2 = this->simplex_index(i, 2);
        double Li1 = a(i1) + b(i1) * x + c(i1) * y + d(i1) * z;
        double Li2 = a(i2) + b(i2) * x + c(i2) * y + d(i2) * z;
        Vector3d bv = basis_coef(this->basis_global_index[i-1]) * this->l(i) *
                        (Li1 * Vector3d(b(i2), c(i2), d(i2)) -
                         Li2 * Vector3d(b(i1), c(i1), d(i1)));
        field += bv;
    }
    return field;
}



