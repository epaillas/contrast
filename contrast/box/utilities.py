import sys
import numpy as np
from scipy.io import FortranFile
import argparse

def ascii_to_unformatted(input_filename, output_filename,
  pos_cols=[0, 1, 2], vel_cols=None, weight_cols=None
):
  # import data
  data = np.genfromtxt(input_filename)
  
  pos = data[:, pos_cols]
  cout = pos # default catalogue with only positions

  if vel_cols is not None:
    vel = data[:, vel_cols]
    cout = np.c_[cout, vel]
  if weight_cols is not None:
    weight = data[:, weight_cols]
    cout = np.c_[cout, weight]

  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()