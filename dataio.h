

#ifndef DATAIO_H_
#define DATAIO_H_

#include "eigen.h"

#include <vector>

struct scalar_field {
    Vector3d posi;
    double value;
};

struct vector_field {
    Vector3d posi;
    Vector3d value;
};

void pos_scalar_field(std::string filename, const std::vector<scalar_field> &dist);
void pos_vector_field(std::string filename, const std::vector<vector_field> &dist);

void export_scalar_field(std::string filename = "scalar_field.pos");
void export_vector_field(std::string filename = "vector_field.pos");

void export_M_matrix();
void export_K_matrix();

std::string integer_to_string(int n);

struct signal {
    
    std::vector<double> ampl;
    double dt;
    
    void push_back(double value) {
        ampl.push_back(value);
    }
    
    void save(std::string filename);

};

#endif /* DATAIO_H_ */
