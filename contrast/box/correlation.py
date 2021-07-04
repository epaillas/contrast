import subprocess
from os import path
import numpy as np


def tpcf_monopole(
    data_filename1, output_filename,
    box_size, dim1_min, dim1_max,
    dim1_nbin, ngrid, data_filename2=None,
    nthreads=1
):
    '''
    Two-point correlation function as a function of the radial
    pairwise separation r.

    Parameters:  data_filename1: str
                 Name of the file containing the data catalogue.

                 output_filename: str
                 Name of the output file that will contain the
                 correlation function.

                 box_size: float
                 Size of the simulation box.

                 dim1_min: float
                 Minimum value for the radial binning.

                 dim1_max: float
                 Maximum value for the radial binning.

                 dim1_nbin: float
                 Number of radial bins.

                 ngrid: int
                 Number of cell divisions to speed up calculations.
                 It has to be a divisor of box_size (e.g. 100 for a
                 1000 Mpc/h box).

                 data_filename2: optional, str
                 Name of the file containing a second data catalogue
                 to perform a cross-correlation function.

                 nthreads: optional, int, defaults to 1
                 Number of threads for parallelization.

    '''

    # if no second data file is provided, assume it
    # is an autocorrelation function.
    if data_filename2 is None:
        data_filename2 = data_filename1

    # check if files exist
    for filename in [data_filename1, data_filename2]:
        if not path.isfile(filename):
            raise FileNotFoundError('{} does not exist.'.format(filename))

    binpath = path.join(path.dirname(__file__),
                        'bin', 'tpcf_monopole.exe')

    cmd = [
        binpath, data_filename1, data_filename2,
        output_filename, str(box_size), str(dim1_min),
        str(dim1_max), str(dim1_nbin), str(ngrid),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')
    subprocess.call(cmd, stdout=log, stderr=log)

    # open output file
    data = np.genfromtxt(output_filename)
    r = data[:, 0]
    corr = data[:, 1]

    return r, corr


def tpcf_monopole_2d(
  data_filename1, output_filename,
  box_size, dim1_min, dim1_max,
  dim1_nbin, ngrid, data_filename2=None,
  nthreads=1
):
    '''
    Projected two-point correlation function as a function of the
    radial pairwise separation r.

    Parameters:  data_filename1: str
                 Name of the file containing the data catalogue.

                 output_filename: str
                 Name of the output file that will contain the
                 correlation function.

                 box_size: float
                 Size of the simulation box.

                 dim1_min: float
                 Minimum value for the radial binning.

                 dim1_max: float
                 Maximum value for the radial binning.

                 dim1_nbin: float
                 Number of radial bins.

                 ngrid: int
                 Number of cell divisions to speed up calculations.
                 It has to be a divisor of box_size (e.g. 100 for a
                 1000 Mpc/h box).

                 data_filename2: optional, str
                 Name of the file containing a second data catalogue
                 to perform a cross-correlation function.

                 nthreads: optional, int, defaults to 1
                 Number of threads for parallelization.

    '''

    # if no second data file is provided, assume it
    # is an autocorrelation function.
    if data_filename2 is None:
        data_filename2 = data_filename1

    # check if files exist
    for filename in [data_filename1, data_filename2]:
        if not path.isfile(filename):
            raise FileNotFoundError('{} does not exist.'.format(filename))

    binpath = path.join(path.dirname(__file__),
                        'bin', 'tpcf_monopole_2d.exe')

    cmd = [
        binpath, data_filename1, data_filename2,
        output_filename, str(box_size), str(dim1_min),
        str(dim1_max), str(dim1_nbin), str(ngrid),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')
    subprocess.call(cmd, stdout=log, stderr=log)

    # open output file
    data = np.genfromtxt(output_filename)
    r = data[:, 0]
    corr = data[:, 1]

    return r, corr


def tpcf_rmu(
    data_filename1, output_filename,
    box_size, dim1_min, dim1_max,
    dim1_nbin, ngrid, data_filename2=None,
    dim2_min=-1, dim2_max=1, dim2_nbin=80,
    nthreads=1
):
    '''
    Two-point correlation function binned in r and mu,
    where r is the pairwise radial separation and mu is
    the cosine of the angle between the pair separation and
    the line of sight.

    Parameters:  data_filename1: str
                 Name of the file containing the data catalogue.

                 output_filename: str
                 Name of the output file that will contain the
                 correlation function.

                 box_size: float
                 Size of the simulation box.

                 dim1_min: float
                 Minimum value for the radial binning.

                 dim1_max: float
                 Maximum value for the radial binning.

                 dim1_nbin: float
                 Number of radial bins.

                 ngrid: int
                 Number of cell divisions to speed up calculations.
                 It has to be a divisor of box_size (e.g. 100 for a
                 1000 Mpc/h box).

                 data_filename2: optional, str
                 Name of the file containing a second data catalogue
                 to perform a cross-correlation function.

                 dim2_min: optional, float, defaults to -1
                 Minimum value for the mu binning.

                 dim2_max: optional, float, defaults to 1
                 Maximum value for the mu binning.

                 dim2_nbin: optional, int, defaults to 80
                 Number of mu bins.

                 nthreads: optional, int, defaults to 1
                 Number of threads for parallelization.

    '''

    # if no second data file is provided, assume it
    # is an autocorrelation function.
    if data_filename2 is None:
        data_filename2 = data_filename1

    # check if files exist
    for filename in [data_filename1, data_filename2]:
        if not path.isfile(filename):
            raise FileNotFoundError('{} does not exist.'.format(filename))

    binpath = path.join(path.dirname(__file__),
                        'bin', 'tpcf_rmu.exe')

    cmd = [
        binpath, data_filename1, data_filename2,
        output_filename, str(box_size), str(dim1_min),
        str(dim1_max), str(dim1_nbin), str(dim2_min),
        str(dim2_max), str(dim2_nbin), str(ngrid),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')
    subprocess.call(cmd, stdout=log, stderr=log)

    # TODO return tpcf array

    return


def mean_radial_velocity_monopole(
    data_filename1, output_filename,
    box_size, dim1_min, dim1_max,
    dim1_nbin, ngrid, data_filename2=None,
    dim2_min=-1, dim2_max=1, dim2_nbin=80,
    nthreads=1
):
    '''
    Mean radial velocity as a function of the radial pairwise
    separation r.

    Parameters:  data_filename1: str
                 Name of the file containing the data catalogue.

                 output_filename: str
                 Name of the output file that will contain the
                 output profile.

                 box_size: float
                 Size of the simulation box.

                 dim1_min: float
                 Minimum value for the radial binning.

                 dim1_max: float
                 Maximum value for the radial binning.

                 dim1_nbin: float
                 Number of radial bins.

                 ngrid: int
                 Number of cell divisions to speed up calculations
                 It has to be a divisor of box_size (e.g. 100 for a
                 1000 Mpc/h box).

                 data_filename2: optional, str
                 Name of the file containing a second data catalogue
                 to perform a cross-correlation function.

                 nthreads: optional, int, defaults to 1
                 Number of threads for parallelization.

    '''

    # if no second data file is provided, assume it
    # is an autocorrelation function.
    if data_filename2 is None:
        data_filename2 = data_filename1

    # check if files exist
    for filename in [data_filename1, data_filename2]:
        if not path.isfile(filename):
            raise FileNotFoundError('{} does not exist.'.format(filename))

    binpath = path.join(path.dirname(__file__),
                        'bin', 'mean_radial_velocity_monopole.exe')

    cmd = [
        binpath, data_filename1, data_filename2,
        output_filename, str(box_size), str(dim1_min),
        str(dim1_max), str(dim1_nbin), str(ngrid),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')
    subprocess.call(cmd, stdout=log, stderr=log)

    # TODO return radial velocity array

    return


def std_los_velocity_rmu(
    data_filename1, output_filename,
    box_size, dim1_min, dim1_max,
    dim1_nbin, ngrid, data_filename2=None,
    dim2_min=-1, dim2_max=1, dim2_nbin=80,
    nthreads=1
):
    '''
    Line-of-sight velocity dispersion inned in r and mu,
    where r is the pairwise radial separation and mu is
    the cosine of the angle between the pair separation and
    the line of sight.

    Parameters:  data_filename1: str
                 Name of the file containing the data catalogue.

                 output_filename: str
                 Name of the output file that will contain the
                 velocity profile

                 box_size: float
                 Size of the simulation box.

                 dim1_min: float
                 Minimum value for the radial binning.

                 dim1_max: float
                 Maximum value for the radial binning.

                 dim1_nbin: float
                 Number of radial bins.

                 ngrid: int
                 Number of cell divisions to speed up calculations.
                 It has to be a divisor of box_size (e.g. 100 for a
                 1000 Mpc/h box).

                 data_filename2: optional, str
                 Name of the file containing a second data catalogue
                 to perform a cross-correlation function.

                 dim2_min: optional, float, defaults to -1
                 Minimum value for the mu binning.

                 dim2_max: optional, float, defaults to 1
                 Maximum value for the mu binning.

                 dim2_nbin: optional, int, defaults to 80
                 Number of mu bins.

                 nthreads: optional, int, defaults to 1
                 Number of threads for parallelization.

    '''

    # if no second data file is provided, assume it
    # is an autocorrelation function.
    if data_filename2 is None:
        data_filename2 = data_filename1

    # check if files exist
    for filename in [data_filename1, data_filename2]:
        if not path.isfile(filename):
            raise FileNotFoundError('{} does not exist.'.format(filename))

    binpath = path.join(path.dirname(__file__),
                        'bin', 'std_los_velocity_rmu.exe')

    cmd = [
        binpath, data_filename1, data_filename2,
        output_filename, str(box_size), str(dim1_min),
        str(dim1_max), str(dim1_nbin), str(dim2_min),
        str(dim2_max), str(dim2_nbin), str(ngrid),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')
    subprocess.call(cmd, stdout=log, stderr=log)

    # TODO return velocity dispersion array

    return
