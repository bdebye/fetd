
#include "dataio.h"
#include "solve.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

const string POS_HEADER = "View \"field\" {";
const string POS_FOOTER = "};";

const char SCALAR_POINT_FORMAT[] = "SP(%f, %f, %f){%e};";
const char VECTOR_POINT_FORMAT[] = "VP(%f, %f, %f){%e, %e, %e};";
char BUFFER[500];

void pos_scalar_field(string filename, const vector<scalar_field> &dist) {
    ofstream ofs(filename);
    ofs << POS_HEADER << endl;
    for(auto iter = dist.begin(); iter != dist.end(); ++iter) {
        sprintf(BUFFER, SCALAR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
                iter->value);
        ofs << BUFFER << endl;
    }
    ofs << POS_FOOTER << endl;
}

void pos_vector_field(string filename, const vector<vector_field> &dist) {
    ofstream ofs(filename);
    ofs << POS_HEADER << endl;
    for(auto iter = dist.begin(); iter != dist.end(); ++iter) {
        sprintf(BUFFER, VECTOR_POINT_FORMAT, iter->posi(0), iter->posi(1), iter->posi(2),
                iter->value(0), iter->value(1), iter->value(2));
        ofs << BUFFER << endl;
    }
    ofs << POS_FOOTER << endl;
}

void signal::save(std::string filename) {
    ofstream ofs(filename);
    std::cout << "Save the received signal..." << std::endl;
    for(int i = 0; i < ampl.size(); i++) {
        ofs << i * dt << "\t" << ampl[i] << endl;
    }
}

int find_max_index(const vector<scalar_field> &dist) {
    double max_v = dist[0].value;
    int max_i = -1;
    for(int i = 0; i < dist.size(); i++) {
        if(dist[i].value > max_v) {
            max_v = dist[i].value;
            max_i = i;
        }
    }
    return max_i;
}

void export_scalar_field(std::string filename) {
    vector<scalar_field> dist;
    dist.reserve(tetrah_element.size());
    //std::cout << tetrah_element[1637].center_field()(2) << std::endl;
    for(int i = 0; i < tetrah_element.size(); i++) {
        scalar_field scalar_f;
        scalar_f.posi = tetrah_element[i].center_coord();
        scalar_f.value = tetrah_element[i].center_field()(2);
        dist.push_back(scalar_f);
    }
    pos_scalar_field(filename, dist);
}

void export_vector_field(std::string filename) {
    vector<vector_field> dist;
    dist.reserve(tetrah_element.size());
    //std::cout << tetrah_element[1637].center_field()(2) << std::endl;
    for(int i = 0; i < tetrah_element.size(); i++) {
        vector_field vector_f;
        vector_f.posi = tetrah_element[i].center_coord();
        vector_f.value = tetrah_element[i].center_field();
        dist.push_back(vector_f);
    }
    pos_vector_field(filename, dist);
}

void export_M_matrix() {
    ofstream ofs("M.txt");
    ofs << MatrixXd(M) << endl;
    ofs.close();
}

void export_K_matrix() {
    ofstream ofs("K.txt");
    ofs << MatrixXd(K) << endl;
    ofs.close();
}

std::string integer_to_string(int n) {
    stringstream ss;
    ss << n;
    return ss.str();
}

