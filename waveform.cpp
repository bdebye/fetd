
#include "waveform.h"

#include <cmath>

double gauss_waveform(double t, double tau) {
    return exp(-pow(t, 2) / pow(tau, 2) / 2);
}

double gauss_differential(double t, double tau) {
    return - (t / pow(tau, 2.0)) * exp(-pow(t, 2.0) / pow(tau, 2.0) / 2.0);
}

double gauss_dd(double t, double tau) {
    return ((pow(t, 2) - pow(tau, 2)) / pow(tau, 4)) * exp(-pow(t, 2) / pow(tau, 2) / 2);
}
