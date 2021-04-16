import sys
import numpy as np
from scipy.io import FortranFile
import argparse

def fits_to_unformatted(
  input_filename, output_filename, omega_m
):
  # define cosmology
  cosmo = Cosmology(omega_m=omega_m)

  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  # convert redshifts to distances
  dist = cosmo.ComovingDistance(cat['Z'])
  x = dist * np.sin(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.sin(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.cos(cat['DEC'] * np.pi / 180)

  # write result to output file
  cout = np.c_[x, y, z]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()

def ascii_to_unformatted(input_filename, output_filename
):
  # import data
  cout = np.genfromtxt(input_filename)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()