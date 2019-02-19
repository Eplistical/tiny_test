#include "misc/ioer.hpp"
#include "misc/vector.hpp"
#include "boris_integrator.hpp"

using namespace std;

const double c = -1.0;
const double m = 1.0;
const double dt = 0.005;
const int Nstep = 10000;

vector<double> cal_w(const vector<double>& r) {
    //return c / m * vector<double> { 0.0, 0.0, 1.0 };
    return c / m * vector<double> { 0.0, 0.0, 1.0 / (r[0]*r[0] + r[1]*r[1]) };
}

vector<double> cal_force(const vector<double>& r) {
    return vector<double> { 0.0, 0.0, 0.0 };
}

int main(int argc, char** argv) {
    vector<double> r { 0.01, 0.01, 0.0 };
    vector<double> p { 1.0, 1.0, 0.0 };
    ioer::tabout("#t", "x", "y", "z", "px", "py", "pz");
    for (int istep(0); istep < Nstep; ++istep) {
        ioer::tabout(istep * dt, r, p);
        boris_integrator(r, p, m, dt, cal_force, cal_w);
    }
        
    return 0;
}
