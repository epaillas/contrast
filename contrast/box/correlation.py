from julia.api import Julia
from os import path
import numpy as np


def tpcf_monopole(
    positions1, rbins, 
    box_size, positions2=None
):
    # import Julia modules
    jl = Julia(compiled_modules=False)
    from julia import Main

    module_path = path.join(path.dirname(__file__),
                        'fastmodules', 'paircounts_r.jl')

    jl.eval(f'include("{module_path}")')

    # if no second positions are provided, assume it
    # is an autocorrelation function.
    if positions2 is None:
        positions2 = positions1

    npos1 = len(positions1)
    npos2 = len(positions2)

    # add this so that Julia can recognize these varaibles
    Main.positions1 = positions1.T
    Main.positions2 = positions2.T
    Main.box_size = box_size
    Main.rbins = rbins

    D1D2 = jl.eval("count_pairs_r(positions1, positions2, box_size, rbins)") 

    mean_density = len(positions2) / (box_size**3)
    R1R2 = np.zeros(len(D1D2))

    for i in range(len(rbins) - 1):
        bin_volume = 4 / 3 * np.pi * (rbins[i+1]**3 - rbins[i]**3)
        R1R2[i] = bin_volume * mean_density * npos1

    delta_r = D1D2 / R1R2 - 1

    return delta_r

def projected_tpcf(
    positions1, rbins, 
    box_size, positions2=None
):
    # import Julia modules
    jl = Julia(compiled_modules=False)
    from julia import Main

    module_path = path.join(path.dirname(__file__),
                        'fastmodules', 'paircounts_r.jl')

    jl.eval(f'include("{module_path}")')

    # if no second positions are provided, assume it
    # is an autocorrelation function.
    if positions2 is None:
        positions2 = positions1

    npos1 = len(positions1)
    npos2 = len(positions2)

    # add this so that Julia can recognize these varaibles
    Main.positions1 = positions1.T
    Main.positions2 = positions2.T
    Main.box_size = box_size
    Main.rbins = rbins

    D1D2 = jl.eval("count_pairs_2d_r(positions1, positions2, box_size, rbins)") 

    mean_density = len(positions2) / (box_size**2)
    R1R2 = np.zeros(len(D1D2))

    for i in range(len(rbins) - 1):
        bin_volume = np.pi * (rbins[i+1]**2 - rbins[i]**2)
        R1R2[i] = bin_volume * mean_density * npos1

    delta_r = D1D2 / R1R2 - 1

    return delta_r
