import sys
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
  Two-point cross-correlation function between
  two sets of points. If only one sample of points
  is provided, then it is equivalent to the two-point
  autocorrelation function.

  Input arguments:

  data_filename1: name of text file containing first
  set of tracers, where the first three columns correspond
  to the positions in cartesian coordinates.

  output_filename: name of text file where the output 
  is going to be stored.

  box_size: size of the simulation box.

  dim1_min: minimum pair distance for the correlation function.

  dim1_max: maximum pair distance for the correlation function.

  dim1_nbin: number of radial bins for the correlation function.

  ngrid: number of cells in which to divide the simulation box 
  to speed up calculation. It has to be a divisor of box_size 
  (e.g. 100 for 1000 Mpc/h box, or 128 for a 1024 Mpc/h box).

  data_filename2 (optional): name of text file containing second 
  set of points for a cross-correlation function.

  nthreads (optional): number of threads to use while running the 
  algorithm.
  '''

  # check if files exist
  if not path.isfile(data_filename1):
    raise FileNotFoundError('{} does not exist.'.format(data_filename1))

  if data_filename2 == None:
    data_filename2 = data_filename1

  binpath = path.join(path.dirname(__file__),
  'bin', 'tpcf.exe')

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
  r = data[:,0]
  corr = data[:,1]

  return r, corr


def tpcf_2d(
  data_filename1, output_filename,
  box_size, dim1_min, dim1_max,
  dim1_nbin, ngrid, data_filename2=None, 
):
  '''
  Projected (2D) two-point cross-correlation function
  between two sets of points. If only one sample of points
  is provided, then it is equivalent to the projected two-point
  autocorrelation function.

  Input arguments:

  data_filename1: name of text file containing first
  set of tracers, where the first three columns correspond
  to the positions in cartesian coordinates.

  output_filename: name of text file where the output 
  is going to be stored.

  box_size: size of the simulation box.

  dim1_min: minimum pair distance for the correlation function.

  dim1_max: maximum pair distance for the correlation function.

  dim1_nbin: number of radial bins for the correlation function.

  ngrid: number of cells in which to divide the simulation box 
  to speed up calculation. It has to be a divisor of box_size 
  (e.g. 100 for 1000 Mpc/h box, or 128 for a 1024 Mpc/h box).

  data_filename2 (optional): name of text file containing second 
  set of points for a cross-correlation function.
  '''

  # check if files exist
  if not path.isfile(data_filename1):
    raise FileNotFoundError('{} does not exist.'.format(data_filename1))

  if data_filename2 == None:
    data_filename2 = data_filename1

  binpath = path.join(path.dirname(__file__),
  'bin', 'tpcf_2d.exe')

  print(binpath)

  cmd = [
    binpath, data_filename1, data_filename2,
    output_filename, str(box_size), str(dim1_min),
    str(dim1_max), str(dim1_nbin), str(ngrid)
  ]

  log_filename = '{}.log'.format(output_filename)
  log = open(log_filename, 'w+')
  subprocess.call(cmd, stdout=log, stderr=log)

  # open output file
  data = np.genfromtxt(output_filename)
  r = data[:,0]
  corr = data[:,1]

  return r, corr

def mean_radial_velocity_monopole(
  data_filename1, output_filename,
  box_size, dim1_min, dim1_max,
  dim1_nbin, ngrid, data_filename2=None, 
  nthreads=1
):
  '''
  Mean radial velocity as a function of the radial pair
  separation.

  Input arguments:

  data_filename1: name of text file containing first
  set of tracers, where the first three columns correspond
  to the positions in cartesian coordinates.

  output_filename: name of text file where the output 
  is going to be stored.

  box_size: size of the simulation box.

  dim1_min: minimum pair distance for the correlation function.

  dim1_max: maximum pair distance for the correlation function.

  dim1_nbin: number of radial bins for the correlation function.

  ngrid: number of cells in which to divide the simulation box 
  to speed up calculation. It has to be a divisor of box_size 
  (e.g. 100 for 1000 Mpc/h box, or 128 for a 1024 Mpc/h box).

  data_filename2 (optional): name of text file containing second 
  set of points for a cross-correlation function.

  nthreads (optional): number of threads to use while running the 
  algorithm.
  '''

  # check if files exist
  if not path.isfile(data_filename1):
    raise FileNotFoundError('{} does not exist.'.format(data_filename1))

  if data_filename2 == None:
    data_filename2 = data_filename1

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

  # open output file
  data = np.genfromtxt(output_filename)
  r = data[:,0]
  corr = data[:,1]

  return r, corr

def tpcf_rmu(
    data_filename1, data_filename2, output_filename,
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
    if data_filename2 == None:
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


def std_los_velocity_rmu(
    data_filename1, data_filename2, output_filename,
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
    if data_filename2 == None:
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
