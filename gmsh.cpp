

#include "gmsh.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

mesh_data mesh;
std::map<int, medium_property> medium_conf;

void clear_mesh() {
    mesh.node_num = 0;
    mesh.triang_num = 0;
    mesh.tetrah_num = 0;
    mesh.segment_num = 0;
    
    mesh.node.clear();
    mesh.segment.clear();
    mesh.triang.clear();
    mesh.tetrah.clear();
}

void read_node(string line) {
    istringstream iss(line);
    Vector3d node;
    int n;
    iss >> n;
    iss >> node(0);
    iss >> node(1);
    iss >> node(2);
    mesh.node.push_back(node);
}

const int SEGMENT_TYPE = 1;
const int TRIANG_TYPE = 2;
const int TETRAH_TYPE = 4;

void read_element(string line) {
    istringstream iss(line);
    int n;
    int code;
    iss >> n;
    iss >> code;
    iss >> n;
    if(code == SEGMENT_TYPE) {
        msh_segment segment;
        iss >> segment.physical_tag;
        iss >> segment.geometrical_tag;
        iss >> segment.node_number_list[0];
        iss >> segment.node_number_list[1];
        sort(segment.node_number_list.begin(), segment.node_number_list.end());
        mesh.segment.push_back(segment);
        mesh.segment_num++;
    }
    else if(code == TRIANG_TYPE) {
        msh_triang triang;
        iss >> triang.physical_tag;
        iss >> triang.geometrical_tag;
        iss >> triang.node_number_list[0];
        iss >> triang.node_number_list[1];
        iss >> triang.node_number_list[2];
        sort(triang.node_number_list.begin(), triang.node_number_list.end());
        mesh.triang.push_back(triang);
        mesh.triang_num++;
    }
    else if(code == TETRAH_TYPE) {
        msh_tetrah tetrah;
        iss >> tetrah.physical_tag;
        iss >> tetrah.geometrical_tag;
        iss >> tetrah.node_number_list[0];
        iss >> tetrah.node_number_list[1];
        iss >> tetrah.node_number_list[2];
        iss >> tetrah.node_number_list[3];
        sort(tetrah.node_number_list.begin(), tetrah.node_number_list.end());
        mesh.tetrah.push_back(tetrah);
        mesh.tetrah_num++;
    }
    
}

void load_gmsh_file(string filename) {
    fstream fs(filename);
    string line;
    istringstream iss;
    int element_num = 0;
    clear_mesh();
    
    if(!fs.good()) {
        cout << "Broken mesh file or file doesn't exist..." << endl;
    }

    mesh.node.push_back(Vector3d());
    for(int i = 0; i < 4; i++)
        getline(fs, line);
    getline(fs, line);
    
    iss.str(line);
    iss >> mesh.node_num;
    mesh.node.reserve(mesh.node_num);
    for(int i = 0; i < mesh.node_num; i++) {
        getline(fs, line);
        read_node(line);
    }
    
    getline(fs, line);
    getline(fs, line);
    getline(fs, line);
    
    iss.clear();
    iss.str(line);
    iss >> element_num;

    for (int i = 0; i < element_num; i++) {
        getline(fs, line);
        read_element(line);
    }
}
