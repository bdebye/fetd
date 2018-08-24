
#include <iostream>
#include <vector>
#

#include "dataio.h"
#include "eigen.h"
#include "element.h"
#include "gmsh.h"
#include "solve.h"
#include "feed.h"

using namespace std;

void print_global_edge_index() {
    for(int e = 0; e < tetrah_element.size(); e++) {
        for(int i = 0; i < 6; i++) {
            cout << tetrah_element[e].basis_global_index[i] << " ";
        }
        cout << endl;
    }
}

void find_center_element() {
    for(int i = 0; i < tetrah_element.size(); i++) {
        if(norm(tetrah_element[i].center_coord() - Vector3d(2, 2, 2)) < 0.23)
            cout << i << endl;
    }
}

int main() {
    load_gmsh_file("fetd_dipole.msh");
    parse_edge_basis();
    generate_elements();
    //generate_boundary_elements();
    
    //for(int i = 0; i < 3; i++) {
     //   cout << boundary_element[0].vertex[i] << endl;
    //}
    
    //cout << boundary_element[0].b(1) << endl;
    //cout << boundary_element[0].b_1_Jin() << endl;
    //assembly();
    //newmark(7000);
    cout << mesh.node_num << endl;
    cout << mesh.segment_num << endl;
    cout << mesh.triang_num << endl;
    cout << mesh.tetrah_num << endl;
}