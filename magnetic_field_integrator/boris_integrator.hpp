#include "misc/vector.hpp"

namespace {

    using std::vector;

    // function type to calculate electric field force 
    using cal_force_type = vector<double>(*)(const vector<double>&);

    // function type to calculate magnetic field
    //  all constant should be included, i.e. F_mag = cross(p, w)
    using cal_w_type = vector<double>(*)(const vector<double>&);

    void boris_integrator(  vector<double>& r, vector<double>& p, 
                            const double m, const double dt,
                            cal_force_type cal_force, cal_w_type cal_w) 
    {
        /*
         * Boris integrator
         *  integerate a particle moving in a magnetic field
         */
        const double half_dt = 0.5 * dt;

        r += half_dt / m * p;
        vector<double> F = cal_force(r);
        p += half_dt * F;

        // p -> R(0.5*dt*omega) * p
        //  R(z) = (I + z) * inv(I - z)
        vector<double> w = cal_w(r) * half_dt;
        const double w0 = w[0], w1 = w[1], w2 = w[2];
        const double w00 = w0*w0, w11 = w1*w1, w22 = w2*w2;
        const double w01 = w0*w1, w02 = w0*w2, w12 = w1*w2;
        const double p0 = p[0], p1 = p[1], p2 = p[2];
        p[0] = (w00 - w11 - w22 + 1.0) * p0 + 2 * (w2 + w01) * p1 + 2 * (w02 - w1) * p2;
        p[1] = 2 * (w01 - w2) * p0 + (w00 + w11 - w22 + 1.0) * p1 + 2 * (w12 + w0) * p2;
        p[2] = 2 * (w1 + w02) * p0 + 2 * (w12 - w0) * p1 + (w00 - w11 + w22 + 1.0) * p2;
        p /= (w00 + w11 + w22 + 1.0);

        p += half_dt * F;
        r += half_dt / m * p;
    }
};
