
#ifndef GMSH_H_
#define GMSH_H_


#include <vector>
#include <string>
#include <array>
#include <map>

#include "eigen.h"


struct msh_segment {
    std::array<int, 2> node_number_list;
    int physical_tag;
    int geometrical_tag;
};

struct msh_triang {
    std::array<int, 3> node_number_list;
    int physical_tag;
    int geometrical_tag;
};

struct msh_tetrah {
    std::array<int, 4> node_number_list;
    int physical_tag;
    int geometrical_tag;
};

struct medium_property {
    double epsilon;
    double sigma;
    double mu;
};


struct mesh_data {
    
    int triang_num;
    int tetrah_num;
    int segment_num;
    int node_num;
    
    std::vector<Vector3d> node;
    std::vector<msh_segment> segment;
    std::vector<msh_triang> triang;
    std::vector<msh_tetrah> tetrah;
    
};

extern mesh_data mesh;
extern std::map<int, medium_property> medium_conf;

void load_gmsh_file(std::string filename);

inline Vector3d get_node(int n) {
    return mesh.node[n];
}

#endif /* defined(GMSH_H_) */
