import glob
import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import simps
from scipy import special


def save_as_unformatted(data, filename):
    '''
    Saves a numpy array as an unformatted
    Fortran 90 file that can be handled by
    this package's numerical routines.

    Parameters:  data: ND array_like
                 Array to be saved.

                 filename: str
                 Name of the output file.
    '''
    data = np.asarray(data)

    nrows, ncols = np.shape(data)
    f = FortranFile(filename, 'w')
    nrows, ncols = np.shape(data)
    f.write_record(nrows)
    f.write_record(ncols)
    f.write_record(data)
    f.close()


def read_array_2d(filename):
    '''
    Read a two dimensional from an ascii
    file.

    Parameters:  filename: str
                 Name of the ascii file containing
                 the array
    '''
    data = np.genfromtxt(filename)
    dim1 = np.unique(data[:, 0])
    dim2 = np.unique(data[:, 1])

    vary_dim2 = False
    if data[0, 0] == data[1, 0]:
        vary_dim2 = True

    result = np.zeros([len(dim1), len(dim2)])
    counter = 0
    if vary_dim2:
        for i in range(len(dim1)):
            for j in range(len(dim2)):
                result[i, j] = data[counter, 2]
                counter += 1
    else:
        for i in range(len(dim2)):
            for j in range(len(dim1)):
                result[j, i] = data[counter, 2]
                counter += 1
    return dim1, dim2, result


def project_to_multipoles(r_c, mu_edges, corr, ells=(0,2,4)):
    if np.ndim(ells) == 0:
        ells = (ells,)
    ells = tuple(ells)
    toret = []
    for ill,ell in enumerate(ells):
        dmu = np.diff(mu_edges, axis=-1)
        poly = special.legendre(ell)(mu_edges)
        legendre = (2*ell + 1) * (poly[1:] + poly[:-1])/2. * dmu
        toret.append(np.sum(corr*legendre, axis=-1)/np.sum(dmu))
    return r_c, toret


def covariance_matrix(data, norm=False):
    '''
    Calculates the covariance matrix of an array
    of measurements.

    Parameters:  data: 2D array-like
                 Array containing the measurements.

                 norm: boolean
                 Set to True for a correlation matrix.
                 False for a covariance matrix.
    '''
    nobs, nbins = np.shape(data)
    mean = np.mean(data, axis=0)
    cov = np.zeros([nbins, nbins])

    for k in range(nobs):
        for i in range(nbins):
            for j in range(nbins):
                cov[i, j] += (data[k, i] - mean[i])*(data[k, j] - mean[j])

    cov /= nobs - 1

    if norm:
        corr = np.zeros_like(cov)
        for i in range(nbins):
            for j in range(nbins):
                corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        return corr
    else:
        return cov


def cross_covariance_matrix(data1, data2, norm=False):
    '''
    Calculates the cross-covariance matrix of two
    sets of measurements.

    Parameters:  data1, data2: 2D array-like
                 Arrays containing the measurements.

                 norm: boolean
                 Set to True for a correlation matrix.
                 False for a covariance matrix.
    '''
    nobs, nbins = np.shape(data1)
    mean1 = np.mean(data1, axis=0)
    mean2 = np.mean(data2, axis=0)
    cov = np.zeros([nbins, nbins])

    for k in range(nobs):
        for i in range(nbins):
            for j in range(nbins):
                cov[i, j] += (data1[k, i] - mean1[i])*(data2[k, j] - mean2[j])

    cov /= nobs - 1

    if norm:
        corr = np.zeros_like(cov)
        for i in range(nbins):
            for j in range(nbins):
                corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        return corr
    else:
        return cov
