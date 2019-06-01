# This python script contains functions and classes used in the
# two synthetic galaxy model notebooks in this directory.

import astropy.io.ascii as asciitable
from scipy import interpolate
import numpy as np


def LoadData(galaxy_datafile, HIrad, ddensdR):
    # Read the galaxy description file which contains rotation velocity and
    # density information in a comma-delimited file where each row contains
    # radius (in kpc), a rotational velocity (in km/s), and a density 
    # (in atoms/cm^3).
    #
    raw_data = asciitable.read(galaxy_datafile)

    # Restructure the data into a 3xN array containing radius, rot. vel.,
    # and gas density in "columns"
    galaxy_data = np.hstack((raw_data['radius'].reshape(-1, 1),
                             raw_data['rot_vel'].reshape(-1, 1),
                             raw_data['density'].reshape(-1, 1)))

    #
    # Interpolate any extension to the rotation curve as flat,
    # but with dropping density out to HI radius limit.
    #
    extrapol_step = (25 - galaxy_data[-1][0])/4
    while (galaxy_data[-1][0] < HIrad):
        new_rad = galaxy_data[-1][0]+extrapol_step
        new_vel = galaxy_data[-1][1]
        if (galaxy_data[-1][2] > ddensdR*extrapol_step):
            new_dens = galaxy_data[-1][2] - ddensdR*extrapol_step
        else:
            new_dens = 0.0
        new_row = np.array([new_rad, new_vel, new_dens])
        galaxy_data = np.vstack((galaxy_data, new_row))

    # Save raw values
    rad_raw = np.copy(galaxy_data[:, 0])
    rotvel_raw = np.copy(galaxy_data[:, 1])
    density_raw = np.copy(galaxy_data[:, 2])

    return (rad_raw, rotvel_raw, density_raw)


def spline_curves(rad, vel, dens, dr):
    # Do a modified spline fit to smooth the rotation curve and gas density
    # data to smooth out data gaps and make sure density and rotational
    # velocity don't go negative

    # Find the spline representation of rotation curve (default is cubic
    # spline)
    rotvel_fit = interpolate.splrep(rad, vel, s=0)
    density_fit = interpolate.splrep(rad, dens, s=0)

    # Fit spline along evenly spaced points (in radius) and restrict rotational
    # velocity and density to be positive (since spline fit is bit wiggly here
    # and caused very small scale 'negative' values at origin for velocity and
    # at high radii for density).
    rad_sp = np.linspace(0, rad[-1], int(rad[-1]/dr))
    rotvel_sp = np.absolute(interpolate.splev(rad_sp, rotvel_fit, der=0).round(1))
    density_sp = np.absolute(interpolate.splev(rad_sp, density_fit, der=0).round(3))

    return(rad_sp, rotvel_sp, density_sp)


def RotCoord(x, y):
    # Converts (x, y) to (r,theta)
    # Can work on entire arrays of x and y
    return (np.sqrt(x*x+y*y), np.arctan2(y, x))


class nf(float):
    # This class allows floating point numbers to be printed as integers.
    # Based on 
    # http://matplotlib.sourceforge.net/examples/pylab_examples/contour_label_demo.html
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
