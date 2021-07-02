import glob
import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps


def mean_from_mocks(input_handle, output_filename):
    '''
    Averages a Contrast measurement over multiple
    mock realizations.

    Parameters: input_handle: str
                Handle containing a wildcard that is used
                to list all realizations of the measurement.

                output_filename: str
                Name of the output file.
    '''
    mock_files = sorted(glob.glob(input_handle))
    data_list = []

    for mock_file in mock_files:
        data = np.genfromtxt(mock_file)
        data_list.append(data)

    data_list = np.asarray(data_list)
    data_mean = np.nanmean(data_list, axis=0)
    data_mean = np.nan_to_num(data_mean)

    np.savetxt(output_filename, data_mean)


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


def get_multipole(ell, s, mu, xi_smu):
    '''
    Calculate the ell multipole moment
    from the redshift-space correlation
    function xi_smu.

    Parameters:  ell: int
                 Multipole moment to calculate.

                 s: 1D array_like
                 Radial bins of the correlation function.

                 mu: 1D array_like
                 Mu bins of the correlation function, where
                 mu is the cosine of the angle with respect
                 to the line of sight.

                 xi_smu: 2D array_like
                 Correlation function binned in s and mu.
    '''
    multipole = np.zeros(xi_smu.shape[0])
    if mu.min() < 0:
        factor = 2
        mumin = -1
    else:
        factor = 1
        mumin = 0
    for i in range(xi_smu.shape[0]):
        mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3, ext=3)
        xaxis = np.linspace(mumin, 1, 1000)
        lmu = eval_legendre(ell, xaxis)
        yaxis = mufunc(xaxis) * (2 * ell + 1) / factor * lmu
        multipole[i] = simps(yaxis, xaxis)
    return s, multipole


def CovarianceMatrix(data, norm=False):
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


def CrossCovarianceMatrix(data1, data2, norm=False):
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
