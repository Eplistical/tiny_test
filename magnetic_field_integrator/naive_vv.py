import numpy as np

c = -1.0
m = 1.0
Nstep = 10000
dt = 0.005

def cal_B(r, p):
    return np.array([0.0, 0.0, 1.0])


def cal_forceB(r, p):
    B = cal_B(r, p)
    return c / m * np.cross(p, B)


def cal_forceE(r, p):
    return np.array([0.0, 0.0, 0.0])


def integrator(r, p, dt):
    F = cal_forceE(r, p) + cal_forceB(r, p)
    p += 0.5 * dt * F
    r += dt * p / m
    F = cal_forceE(r, p) + cal_forceB(r, p)
    p += 0.5 * dt * F
    return r, p


def run():
    r = np.array([1.0, 0.0, 0.0])
    p = np.array([0.0, 1.0, 0.0])

    print('#%7s%8s%8s%8s%8s%8s%8s' % ('t', 'x', 'y', 'z', 'px', 'py', 'pz'))
    for istep in range(Nstep):
        print('%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f' % (istep * dt, r[0], r[1], r[2], p[0], p[1], p[2]))
        r, p = integrator(r, p, dt)


if __name__ == '__main__':
    run()


# END

