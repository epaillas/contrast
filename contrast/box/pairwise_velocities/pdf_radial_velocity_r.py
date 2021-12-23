from julia.api import Julia
from os import path, environ


def pdf_radial_velocity_r(
    positions1, velocities1, rbins, 
    box_size, positions2=None, velocities2=None,
    nthreads=1
):
    environ['JULIA_NUM_THREADS'] = f'{nthreads}'
    jl = Julia(compiled_modules=False)
    from julia import Main

    module_path = path.join(path.dirname(__file__),
                        'fastmodules', 'pairwise_velocities.jl')

    jl.eval(f'include("{module_path}")')

    # if no second positions are provided, assume it
    # is an autocorrelation function.
    if positions2 is None:
        positions2 = positions1

    if velocities2 is None:
        velocities2 = velocities1

    # add this so that Julia can recognize these varaibles
    Main.positions1 = positions1.T
    Main.positions2 = positions2.T
    Main.velocities1 = velocities1.T
    Main.velocities2 = velocities2.T
    Main.box_size = box_size
    Main.rbins = rbins

    D1D2, V1V2 = jl.eval("pdf_radial_velocity_r(positions1, positions2, velocities1, velocities2, box_size, rbins)") 

    v_r = V1V2 / D1D2

    return v_r

