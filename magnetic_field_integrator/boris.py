import numpy as np
from scipy import linalg as la

c = -1.0
m = 1.0
Nstep = 10000
dt = 0.005


def cal_B(r, p):
    return np.array([0.0, 0.0, 1.0])
    #return np.array([0.0, 0.0, 1.0 / r[0]**2])


def cal_forceB(r, p):
    B = cal_B(r, p)
    return c * np.cross(p, B)


def cal_forceE(r, p):
    return np.array([0.0, 0.0, 0.0])


def cal_R(z):
    I = np.eye(3)
    return np.dot((I + z), la.inv(I - z))


def integrator(r, p, dt):
    r += 0.5*dt * p / m

    F = cal_forceE(r, p)
    p += 0.5*dt * F

    B = cal_B(r, p)
    omega_vec = c / m * B 
    omega_mat = np.array([
                    [               0,  omega_vec[2], -omega_vec[1]], 
                    [   -omega_vec[2],             0,  omega_vec[0]], 
                    [    omega_vec[1], -omega_vec[0],             0]
                    ])
    R = cal_R(0.5*dt * omega_mat)

    p = np.dot(R, p)
    #p += 1.0*dt * c / m * np.cross(p, B) 
    
    p += 0.5*dt * F

    r += 0.5*dt * p / m
    return r, p


def run():
    r = np.array([1.0, 0.0, 0.0])
    p = np.array([0.0, 1.0, 0.0])

    print('#%11s%12s%12s%12s%12s%12s%12s' % ('t', 'x', 'y', 'z', 'px', 'py', 'pz'))
    for istep in range(Nstep):
        I = p[1] + 1.0 / r[0]
        #print('%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f' % (istep * dt, r[0], r[1], r[2], p[0], p[1], p[2]))
        print('%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f%12.2f' % (istep * dt, r[0], r[1], r[2], p[0], p[1], p[2], I))
        r, p = integrator(r, p, dt)


if __name__ == '__main__':
    run()


# END

