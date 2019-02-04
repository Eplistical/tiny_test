#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt

'''
evaluate integral I = \int_0^1 exp(-x)dx with Monte Carlo

compare uniform sampling and importance sampling

-- Assume one can only sample linear functions
'''


def uniform_samp(integrand, Nstep):
    ''' 
    sample function f(x) = 1
    '''
    rst = 0.0

    for i in range(Nstep):
        x = np.random.uniform()
        rst += integrand(x)

    return rst / Nstep



def importance_samp(integrand, Nstep, alpha):
    '''
    sample function f(x) = (1-alpha*x) / (1-alpha/2)
    integrand -> integrand / f(x)
    '''
    rst = 0.0

    for i in range(Nstep):
        # sample f(x), P(x) = int_0^x f(z)dz
        Px = np.random.uniform()
        x = (1.0 - (1.0 - alpha*(2.0-alpha)*Px)**0.5) / alpha
        # cal f(x)
        fx = (1 - alpha*x) / (1-alpha/2)
        rst += integrand(x) / fx

    return rst / Nstep


def run():
    integrand = lambda x : np.exp(-x)
    Nstep = 100000
    Ntrial = 30
    print('# MC calculation for int_0^1 exp(-x)dx: ')
    print('# Ntrial = %d, Nstep per trial = %d' % (Ntrial, Nstep))

    print('#%16s%16s%16s' % ('alpha', 'mean', 'std'))
    print('#%16s%16.8f%16.8f # exact' % ('', 1-np.exp(-1), 0))

    uniform_rst = np.zeros(Ntrial)
    importance_rst = np.zeros(Ntrial)

    for i in range(Ntrial):
        uniform_rst[i] = uniform_samp(integrand, Nstep)
    print(' %16.2f%16.8f%16.8f # uniform' % (0.0, np.mean(uniform_rst), np.std(uniform_rst)))

    for alpha in np.linspace(0.01, 1, 100):
        for i in range(Ntrial):
            importance_rst[i] = importance_samp(integrand, Nstep, alpha)
        print(' %16.2f%16.8f%16.8f # importance' % (alpha, np.mean(importance_rst), np.std(importance_rst)))


if __name__ == '__main__':
    run()
