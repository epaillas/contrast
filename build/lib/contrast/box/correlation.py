import sys
import subprocess
from os import path
import numpy as np

def tpcf(
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

def mean_radial_velocity_vs_r(
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
  'bin', 'mean_radial_velocity_vs_r.exe')

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