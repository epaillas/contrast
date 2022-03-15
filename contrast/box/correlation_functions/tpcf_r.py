from os import path, environ
import numpy as np
from julia.api import Julia


def tpcf_r(
    positions1, rbins, boxsize, 
    positions2=None, weights1=None, weights2=None,
    nthreads=1, return_paircounts=False
):


    if not hasattr(boxsize, "__len__"):
        boxsize = np.array([boxsize, boxsize, boxsize])
    else:
        boxsize = np.asarray(boxsize)

    environ['JULIA_NUM_THREADS'] = f'{nthreads}'
    jl = Julia(compiled_modules=False)
    from julia import Main

    module_path = path.join(path.dirname(__file__),
                        'fastmodules', 'paircounts.jl')

    jl.eval(f'include("{module_path}")')

    if weights1 is None:
        weights1 = np.ones(len(positions1))

    # if no second positions are provided, assume it
    # is an autocorrelation function.
    if positions2 is None:
        positions2 = positions1
        weights2 = weights1
    elif weights2 is None:
        weights2 = np.ones(len(positions2))

    # add this so that Julia can recognize these varaibles
    Main.positions1 = positions1.T
    Main.positions2 = positions2.T
    Main.weights1 = weights1.T
    Main.weights2 = weights2.T
    Main.boxsize = boxsize
    Main.rbins = rbins

    D1D2 = jl.eval("count_pairs_r(positions1, positions2, weights1, weights2, boxsize, rbins)") 

    mean_density = np.sum(weights2) / (boxsize[0]*boxsize[1]*boxsize[2])
    D1R2 = np.zeros(len(D1D2))

    for i in range(len(rbins) - 1):
        bin_volume = 4 / 3 * np.pi * (rbins[i+1]**3 - rbins[i]**3)
        D1R2[i] = bin_volume * mean_density * np.sum(weights1)

    delta_r = D1D2 / D1R2 - 1

    if return_paircounts:
        return delta_r, D1D2, D1R2 
    else:
        return delta_r

