
#include "boundary.h"
#include "solve.h"

const int BDY_GROUP_TAG = 288;

//std::vector<bdy_element> bdy_element;

bdy_element bdy_element::creat_boundary_element(int n) {
    bdy_element bdy_elem;
    bdy_elem.bdy_element_index = n;
    msh_triang trian = mesh.triang[n];
    
    bdy_elem.basis_global_index[0] = find_edge(trian.node_number_list[0], trian.node_number_list[1]);
    bdy_elem.basis_global_index[1] = find_edge(trian.node_number_list[0], trian.node_number_list[2]);
    bdy_elem.basis_global_index[2] = find_edge(trian.node_number_list[1], trian.node_number_list[2]);
    
    bdy_elem.epsilon = 8.854e-012;
    bdy_elem.sigma = 0;
    bdy_elem.mu = 1.257e-006;
    
    Vector3d p1 = get_node(trian.node_number_list[0]);
    Vector3d p2 = get_node(trian.node_number_list[1]);
    Vector3d p3 = get_node(trian.node_number_list[2]);
    
    bdy_elem.vertex[0](0) = 0;
    bdy_elem.vertex[0](1) = 0;
    bdy_elem.vertex[1](0) = norm(p2 - p1);
    bdy_elem.vertex[1](1) = 0;
    bdy_elem.vertex[2](0) = dot(p3 - p1, p2 - p1) / norm(p2 - p1);
    double cosphi1 = bdy_elem.vertex[2](0) / norm(p1 - p3);
    double sinphi1 = sqrt(1.0 - pow(cosphi1, 2));
    bdy_elem.vertex[2](1) = sinphi1 * norm(p3 - p1);
    
    Matrix3d coord_matrix;
    for(int j = 0; j < 3; j++) {
        coord_matrix(0, j) = bdy_elem.vertex[j](0);
        coord_matrix(1, j) = bdy_elem.vertex[j](1);
        coord_matrix(2, j) = 1;
    }
    bdy_elem.area = fabs(coord_matrix.determinant()) / 2.0;
    bdy_elem.coef_matrix = coord_matrix.inverse();
    
    return bdy_elem;
}

double bdy_element::l(int e) {
    int node_a = edge_basis[this->basis_global_index[e-1]][0];
    int node_b = edge_basis[this->basis_global_index[e-1]][1];
    return norm(get_node(node_a) - get_node(node_b));
}

int bdy_element::simplex_index(int e, int i) {
    const static int index[3][2] = {{1, 2}, {1, 3}, {2, 3}};
    return index[e-1][i-1];
}

double bdy_element::F(int i, int j) {
    const static Matrix3d M = (Matrix3d::Ones() + Matrix3d::Identity()) / 6.0;
    int i1 = simplex_index(i, 1);
    int i2 = simplex_index(i, 2);
    int j1 = simplex_index(j, 1);
    int j2 = simplex_index(j, 2);
    return area * l(i) * l(j)
      * (p(i2, j2) * M(i1-1, j1-1)
       - p(i2, j1) * M(i1-1, j2-1)
       - p(i1, j2) * M(i2-1, j1-1)
       + p(i1, j1) * M(i2-1, j2-1));
}

double bdy_element::W(int i, int j) {
    return F(i, j) * Y();
}
