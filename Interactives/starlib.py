#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
StarLib

This is a set of routines we have written to help in the generation of
stars for the astronomical interactives project.

Functions included
------------------
Rad_calc(mass):
    Returns the radius of a star of a given mass.
Temp_calc(mass):
    Returns the temperature of a star of a given mass.
ConfigStar(mass):
    Returns the radius, temperature and hexcolor of a star of a given mass.
StarMesh(temp, rad, scale, pos):
    Returns a pythreejs Mesh object corresponding to a star.
OldStarMesh(temp, rad, scale, pos):
    Returns a pythreejs Mesh object corresponding to a star using old texture
    approach.
xyplane(max_dist, grid_space):
    Returns a pythreejs Mesh and SurfaceGrid corresponding to the xy plane.
origin_marker(size, axis_rad, axis_color):
    Returns 3 pythreejs Mesh objects marking the origin.
axes(max_dist, axis_rad, axis_color):
    Returns 3 pythreejs Mesh objects representing the x, y, and z axes..
axes_line(max_dist):
    Returns a pythreejs Line object representing the XYZ axes.
OrbitalInfo(mass1, mass2, a, e, phi, N):
    Returns orbital period, perihelion, aphelion, and Pandas dataframe time
    series information of orbits.
RadVelInfo(orbit_info, incl):
    Returns radial velocity Pandas dataframe time series information
    corresponding to a given orbit_info Pandas dataframe and inclination angle.
LightCurveInfo(orbit_info, incl, rad1, rad2, temp1, temp2, Na, Ntheta)
    Returns a light curve Pandas dataframe time series information
    corresponding to given orbit_info, inclination, and stellar parameters.

Classes Declared
----------------
BinaryStarModel:
    An object to handle the model of a binary star system, computing and
    storing orbital information, radial velocity, and light curve information.
BinaryStarViewer:
    An object to handle creation and updating of a pythreejs 3D Render view
    of a BinaryStarModel.

Created on Sun Jun 17 10:43:24 2018

@author: Juan Cabanela
(although many routines were originally written by Sam Holen and
Andrew Louwagie Gordon and just consolidated here)
"""

import numpy as np
import pandas as pd
import pythreejs as p3j
import tempNcolor as tc
import traitlets

# Define common units (in SI units)
hr = 3600             # An hour in seconds
day = 24*hr           # A day in seconds
yr = 3.15581450e7     # A year in seconds
AU = 1.4959787066e11  # An astronomical unit in meters
pc = 206264.806 * AU  # A parsec in meters
deg2rad = np.pi/180   # Multiplier to convert degress to radians
rad2deg = 180/np.pi   # Multiplier to convert radians to degrees


# Define physical constants as globals (in SI units)
G = 6.673e-11        # Universal gravitational constant (in N*m^2/kg^2)
c = 2.99792458e08    # Speed of light (in m/s)
sigma = 5.670367e-8  # Stefan-Boltzmann constant (in W/(m^2*K^4))
M_Sun = 1.9891e30    # Mass of Sun (in kg)
S_Sun = 1.365e3      # Solar constant (in W/m^2)
L_Sun = 4*np.pi*AU*AU*S_Sun  # Solar luminosity (in W)
R_Sun = 6.95508e8    # Solar radius (in m)
# Effective temperature of Sun (in K))
Te_Sun = pow(L_Sun/(4*np.pi*R_Sun*R_Sun*sigma), 0.25)


def Rad_calc(mass=1):
    """
    Estimates radius of a star given the mass.

    Uses approximate fit from ZAMS line from Eker et al (2015) with a small
    tweak to ensure a solar mass star has a radius of 1 solar radii.


    Parameters
    ----------
    mass : float
            Mass of star in solar masses (default 1.0)

    Returns
    -------
    radius : float
              radius of star in units of solar radius.
    """
    return 10 ** (0.0757*(np.log10(mass))**4
                  - 0.1348*(np.log10(mass))**3
                  - 0.1355*(np.log10(mass))**2
                  + 0.8546*np.log10(mass))


def Temp_calc(mass=1):
    """
    Estimates the approximate temperature of a
    star given its mass. Uses a fit aquired from Eker et al (2015)
    with an adjustment to ensure Sun has a temperature of 5777K.

    Parameters
    ----------
    mass : float
            Mass of star in solar masses (default 1.0)

    Returns
    -------
    temperature : float
                   temperature of star in Kelvin
    """
    return 10 ** (-0.0044*np.log10(mass)**4
                  - 0.1600*np.log10(mass)**3
                  + 0.2300*np.log10(mass)**2
                  + 0.6084*np.log10(mass) + 3.7617)


def ConfigStar(mass=1.0):
    """
    Determines the radius (in solar radii), temperature (in K), and
    hexcolor of a star of given mass, assuming it is a main sequence star.

    Parameters
    ----------
    mass : float
            Mass of star in solar masses (default 1.0)

    Returns
    -------
    radius : float
              radius of star in units of solar radius.
    temp : float
              temperature of star in Kelvin
    hexcolor: string
              hexcolor corresponding to stellar temperature.
    """

    # Determine approximate radius in solar radii for both stars.
    radius = Rad_calc(mass)

    # Determines the approximate temperature of each star.
    temp = Temp_calc(mass)

    # Use scalar temperature to estimate hexcolor appropriate to each star
    hexcolor = tc.rgb2hex(tc.temp2rgb(temp))[0]

    return (radius, temp, hexcolor)


def OldStarMesh(temp=Te_Sun, rad=1,  scale=(1, 1, 1), pos=[0, 0, 0]):
    """
    This function creates a pythreejs object that represents a star using
    a texture based on public domain STEREO Heliographic map made with
    images taken of the Sun on Dec. 30, 2011.  Image downloaded from
       https://stereo.gsfc.nasa.gov/360blog/
    and had its brightness and resolution rescaled.

    Parameters
    ----------
    temp : float
            temperature of star in Kelvin (default 5777)
    rad : float
            radius of the star in system units (default 1)
    scale : tuple
             pythreejs scale in each dimension as tuple (default (1, 1, 1) )
    pos : list
           three-dimensional position as list (default [0, 0, 0] )

    Returns
    -------
    star : pythreejs.Mesh
              a spherical pythreejs Mesh object representing a star
    """

    # Check is position is a list
    if isinstance(pos, list):
        # Check if this is a list of 3 items
        if (len(pos) != 3):
            raise TypeError('pos passed to StarMesh must be list of 3 numbers')
        # Check that all the items in the list are numbers
        for this_pos in pos:
            try:
                i = float(this_pos)
            except ValueError:
                raise TypeError('ValueError: pos has non-numerical list item.')
    else:
        raise TypeError('pos passed to StarMesh must be list of 3 numbers')

    # Check is scale is a tuple
    if isinstance(scale, tuple):
        if (len(scale) != 3):
            raise TypeError('scale must be tuple of 3 numbers')
    else:
        raise TypeError('scale must be a tuple')

    # Define color and texture of star surface
    hexcolor = tc.rgb2hex(tc.temp2rgb(float(temp)))[0]
    StarTexture = p3j.ImageTexture(imageUri='images/sun_surface.jpg')

    # Create sphere using MeshBasicMaterial (which is unaffected by lighting)
    StarSurface = p3j.MeshBasicMaterial(color=hexcolor, map=StarTexture)

    StarGeom = p3j.SphereBufferGeometry(radius=rad, widthSegments=32,
                                        heightSegments=16)
    return p3j.Mesh(geometry=StarGeom, material=StarSurface,
                    position=pos, scale=scale)


def StarMesh(temp=Te_Sun, rad=1,  scale=(1, 1, 1), pos=[0, 0, 0]):
    """
    This function creates a pythreejs object that represents a star using
    a set of nested spheres that are partly transparent to simulate limb
    darkening.

    Parameters
    ----------
    temp : float
            temperature of star in Kelvin (default 5777)
    rad : float
            radius of the star in system units (default 1)
    scale : tuple
             pythreejs scale in each dimension as tuple (default (1, 1, 1) )
    pos : list
           three-dimensional position as list (default [0, 0, 0] )

    Returns
    -------
    star : pythreejs.Mesh
              a spherical pythreejs Mesh object representing a star
    """

    # Check is position is a list
    if isinstance(pos, list):
        # Check if this is a list of 3 items
        if (len(pos) != 3):
            raise TypeError('pos passed to StarMesh must be list of 3 numbers')
        # Check that all the items in the list are numbers
        for this_pos in pos:
            try:
                i = float(this_pos)
            except ValueError:
                raise TypeError('ValueError: pos has non-numerical list item.')
    else:
        raise TypeError('pos must be list of 3 numbers')

    # Check is scale is a tuple
    if isinstance(scale, tuple):
        if (len(scale) != 3):
            raise TypeError('scale must be tuple of 3 numbers')
    else:
        raise TypeError('scale must be a tuple')

    # Define color of star surface
    hexcolor = tc.rgb2hex(tc.temp2rgb(float(temp)))[0]

    # Number of transparent layers to use to build star
    layers = 20

    # Radial scaling
    drad = 0.97

    # Starting resolution
    N_azi,  N_pol = 32, 16

    # Layer opacity
    tau = 0.4

    # Radii to use
    radii = rad*drad**np.array(range(layers-1, -1, -1))

    # 3D object to represent star
    star = p3j.Object3D()

    # Build the object from inside out
    for i in range(layers):
        # Tweak number of vertices in sphere up for the outer surface sphere
        if (i > (layers-2)):
            N_azi *= 2
            N_pol *= 2

        geom = p3j.SphereBufferGeometry(radius=radii[i],
                                        widthSegments=N_azi,
                                        heightSegments=N_pol,
                                        renderOrder=-i)

        material = p3j.MeshBasicMaterial(color=hexcolor,
                                         transparent=True,
                                         opacity=tau)
        star.add(p3j.Mesh(geom, material))

    # Set the position and scale of this mesh
    star.position = pos
    star.scale = scale

    return star


def StarMeshColor(star, color):
    """
    This function allows you to change the color of a StarMesh in
    one call (since it is composed of multiple sphere objects, it
    can't be done in a single call).

    Parameters
    ----------
    star : pythreejs.Mesh
              a spherical pythreejs Mesh object representing a star
    color : color
             Any pythreejs color
    """

    for i in range(len(star.children)):
        star.children[i].material.color = color


def xyplane(max_dist, grid_space):
    """
    Generates and returns two pythreejs items: a mesh of a flat surface and a
    SurfaceGrid object, both representing the xy plane.

    NOTE: max_dist will be rounded UP to next grid_space position (e.g. if
       yoru submit a max_width of 38 but a grid_space of 5, the max_width will
       be silently rounded up to 40 to allow an integer number of grids).

    Keyword arguments:
    max_dist (float): Maximum extent of grid from origin in each dimension
    grid_space (float): The grid spacing in system units.


    Parameters
    ----------
    max_dist : float
                maximum extent of grid from origin in each dimension
    grid_space : float
            the grid spacing in system units.

    Returns
    -------
    surface : pythreejs.Mesh
                a  pythreejs Mesh object representing the xy plane.
    surf_grid : pythreejs.SurfaceGrid
                a  pythreejs SurfaceGrid object for the grid on the xy plane.
    """

    # Determine the gridsize to use, adding one additional step if necessary
    x_steps = int(np.ceil(max_dist/grid_space))
    xmax = x_steps * grid_space

    # Number of vertices (subtract one for number of boxes shown)
    nx, ny = (2 * x_steps + 1, 2 * x_steps + 1)
    x = np.linspace(-xmax, xmax, nx)
    y = np.linspace(-xmax, xmax, ny)
    xx, yy = np.meshgrid(x, y)
    z = 0*xx + 0*yy

    # Generate the 3D surface and surface grid to return
    surf_g = p3j.SurfaceGeometry(z=list(z[::-1].flat),
                                 width=2*xmax,
                                 height=2*xmax,
                                 width_segments=nx-1,
                                 height_segments=ny-1)
    # Color '#2F4F4F' also known as 'darkslategrey'
    surface_material = p3j.MeshBasicMaterial(color='#2F4F4F',
                                             transparent=True, opacity=0.3)
    surf = p3j.Mesh(geometry=surf_g, material=surface_material)
    # Color '#808080' also known as 'grey'
    grid_material = p3j.LineBasicMaterial(color='#808080')

    # To avoid overlap, lift grid slightly above the plane scaling
    # by size of grid
    surfgrid = p3j.SurfaceGrid(geometry=surf_g, material=grid_material,
                               position=[0, 0, 5e-3*max_dist])

    return surf, surfgrid


def origin_marker(size, axis_rad=0.25, axis_color='#FFFF00'):
    """
    Generate X, Y, Z cylinders of length size in the form of a pythreejs
    Line object.

    Parameters
    ----------
    size : float
                maximum extent of grid from origin in each dimension
    axis_rad : float
                radius of cylinder representing each axis (default: 0.25)
    axis_color : color
                color the axes are drawn in (default: '#FFFF00' aka yellow)

    Returns
    -------
    Xaxis, Yaxis, Zaxis : pythreejs.Mesh*3
            Three pythreejs Mesh objects representing the x, y, and z axes.
    """
    Xaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=size,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[0, 0, 0])
    Xaxis.rotateZ(np.pi/2)

    Yaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=size,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[0, 0, 0])

    Zaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=size,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[0, 0, 0])
    Zaxis.rotateX(np.pi/2)

    return Xaxis, Yaxis, Zaxis


def axes(max_dist, axis_rad=0.25, axis_color='yellow'):
    """
    Generate cylinders of length max_width along the positive X, Y, and Z
    axes in the form of a pythreejs Mesh objects.

    Parameters
    ----------
    max_dist : float
                maximum (positive) extent of axis from origin in each dimension
    axis_rad : float
                radius of cylinder representing each axis (default: 0.25)
    axis_color : color
                color the axes are drawn in (default: 'yellow')

    Returns
    -------
    Xaxis, Yaxis, Zaxis : pythreejs.Mesh*3
            Three pythreejs Mesh objects representing the x, y, and z axes.
    """
    Xaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=max_dist,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[max_dist/2, 0, 0])
    Xaxis.rotateZ(np.pi/2)

    Yaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=max_dist,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[0, max_dist/2, 0])

    Zaxis = p3j.Mesh(geometry=p3j.CylinderBufferGeometry(radiusTop=axis_rad,
                                                         radiusBottom=axis_rad,
                                                         height=max_dist,
                                                         radiusSegments=12,
                                                         heightSegments=1,
                                                         openEnded=False,
                                                         thetaStart=0,
                                                         thetaLength=2*np.pi),
                     material=p3j.MeshBasicMaterial(color=axis_color),
                     position=[0, 0, max_dist/2])
    Zaxis.rotateX(np.pi/2)

    return Xaxis, Yaxis, Zaxis


def axes_lines(max_dist):
    """
    Generate X, Y, Z axes of length max_width in the form of a pythreejs
    Line object.

    Parameters
    ----------
    max_dist : float
                maximum (positive) extent of axis from origin in each dimension

    Returns
    -------
    axes : pythreejs.Line
            a pythreejs Line object representing the xyz axes.
    """
    axes_geom = p3j.Geometry(vertices=[[0, 0, 0], [max_dist, 0, 0],
                                       [0, 0, 0], [0, max_dist, 0],
                                       [0, 0, 0], [0, 0, max_dist]],
                             colors=['white', 'white', 'white',
                                     'white', 'white', 'white'])

    foo = p3j.Line(geometry=axes_geom,
                   material=p3j.LineBasicMaterial(linewidth=1,
                                                  vertexColors='VertexColors'))
    return foo


class BinaryStarModel(traitlets.HasTraits):
    """
    A model of a binary star system, it has the following attributes

    Attributes
    ----------
    mass1 : float
             mass of star 1 in solar masses (default: 1)
    mass2 : float
             mass of star 2 in solar masses (default: 1)
    a : float
         semi-major axis of binary system reduced mass in AU (default: 1)
    e : float
         eccentricity of binary system reduced mass (default: 0)
    phi : float
           angle between semimajor axis and x axis in degrees (default: 0)
    rad1 : float
          radius of star 1 in solar radii (default: main seq. value for mass1)
    rad2 : float
          radius of star 2 in solar radii (default: main seq. value for mass2)
    temp1 : float
          temperature of star 1 in K (default: main sequence value for mass1)
    temp2 : float
          temperature of star 2 in K (default: main sequence value for mass2)
    rv_sys : float
        radial velocity of center of mass of system in km/s (default: 0)
    incl : float
         angle of orbital plane to sky plane in degrees (default: 0)
    P : float
         orbital period in days
    ap : float
          Periastron distance in AU
    aa : float
          Apastron distance in AU
    maxrad : float
          maximum distance of either star from center of mass in AU
    orbit_info : Pandas dataframe
                 time series positions (in AU) and velocities (in km/s) of
                 both stars.
    collision : Boolean
                 indicates if stars collide in this model
    radvel_info : Pandas dataframe
                 time series radial velocities (in km/s) of both stars.
    lc_info : Pandas dataframe
                 time series percentage of total system flux visible.
    N : int
         number of time steps to use in single period. (default: 1000)
    Na : int
          number of annuli used in computing eclipse fraction (default: 100)
    Ntheta : int
          angular steps used in computing eclipse fraction (default: 360)
    up_to_date : Boolean
          Property to indicate if currently stored orbital, radial velocity,
          or light curve models are up to date with current system parameters.
    continuous_update : Boolean
          Determines if object will continuously update initialized orbital,
          radial velocity, and light curve models when attribute changes affect
          them (default : True)
    a_in_AU : Boolean
          Determines if we assume semimajor axis in AU (versus Solar Radii)

    Methods
    -------
    force_update():
        Forces update of the orbital, radial velocity, and light curve models.
    set_orbital_info():
        Assigns values for various parameters based on orbital model.
    set_radvel_info():
        Fills the radial velocity Pandas dataframe, radvel_info.
    set_lc_info():
        Fills the light curve Pandas dataframe, lc_info.
    wipe_radvel_info():
        Set radvel_info to None and turn off computation of radial velocities.
    wipe_lc_info():
       S et lc_info to None and turn off computation of light curve.
    """

    # Define your traits here, as class attributes
    mass1 = traitlets.Float(allow_none=True)
    mass2 = traitlets.Float(allow_none=True)
    a = traitlets.Float(allow_none=True)
    e = traitlets.Float(allow_none=True)
    phi = traitlets.Float(allow_none=True)
    N = traitlets.Int(allow_none=True)
    incl = traitlets.Float(allow_none=True)
    rv_sys = traitlets.Float(allow_none=True)
    rad1 = traitlets.Float(allow_none=True)
    rad2 = traitlets.Float(allow_none=True)
    temp1 = traitlets.Float(allow_none=True)
    temp2 = traitlets.Float(allow_none=True)
    Na = traitlets.Int(allow_none=True)
    Ntheta = traitlets.Int(allow_none=True)
    mdl_counter = traitlets.Integer(allow_none=True)

    def __init__(self, mass1=1, mass2=1, a=1, e=0, phi=0, rv_sys=0,
                 incl=0, rad1=None, rad2=None, temp1=None, temp2=None,
                 N=1000, Na=100, Ntheta=360, rv_init=False, lc_init=False,
                 a_in_AU=True):
        # initialize the class by calling the superclass
        super(BinaryStarModel, self).__init__()

        # Not enough information set to define orbit, rv, or light curve
        self.up_to_date = False
        self.continuous_update = False
        self.mdl_counter = 0  # Total number of models run

        # Initialize radial velocity and light curves info as blank
        self.radvel_info = None
        self.lc_info = None

        # Process variables without setters/getters
        self.mass1 = float(mass1)
        self.mass2 = float(mass2)
        self.a_in_AU = a_in_AU
        self.a = float(a)
        self.rv_sys = float(rv_sys)
        self.N = int(N)
        self.Na = int(Na)
        self.Ntheta = int(Ntheta)

        # Process variables with setters/getters
        self.e = e
        self.phi = phi
        self.incl = incl

        # Initialize radius and temperature values to main sequence
        # values if not set
        if (rad1 is None):
            self.rad1 = Rad_calc(self.mass1)
        else:
            self.rad1 = float(rad1)
        if (rad2 is None):
            self.rad2 = Rad_calc(self.mass2)
        else:
            self.rad2 = float(rad2)
        if (temp1 is None):
            self.temp1 = Temp_calc(self.mass1)
        else:
            self.temp1 = float(temp1)
        if (temp2 is None):
            self.temp2 = Temp_calc(self.mass2)
        else:
            self.temp2 = float(temp2)

        # Call orbital information once quantities initialized
        self.set_orbital_info()

        # Initialize radial velocity data if requested
        if (rv_init):
            self.set_radvel_info()

        # Initialize light curve data if requested
        if (lc_init):
            self.set_lc_info()

        # Indicate up to date and turn on continuous updating
        self.up_to_date = True
        self.continuous_update = True

    ##
    # Define Methods that watch for traitlet definitions and changes
    ##

    @traitlets.validate('a')
    def _check_a(self, proposal):
        # Define semimajor axis as property that gets checked and internally
        # stored in AU even if input for solar radii
        test_a = proposal['value']
        if (self.a_in_AU):
            self.a_AU = test_a
        else:  # input semimajor in solar radii, convert to AU
            self.a_AU = test_a * R_Sun/AU
        return test_a

    @traitlets.validate('e')
    def _check_e(self, proposal):
        test_e = proposal['value']
        # Define eccentricity as 0 through (but not including) 1
        if ((test_e < 0) or (test_e >= 1)):
            raise traitlets.TraitError('e must be >=0 and <1')
        return test_e

    @traitlets.validate('phi')
    def _check_phi(self, proposal):
        # Restrict phi (major axis longitude)
        test_phi = proposal['value']
        if ((test_phi < 0) or (test_phi >= 360)):
            raise traitlets.TraitError('phi must be in range 0 to 360')
        return test_phi

    @traitlets.validate('incl')
    def _check_incl(self, proposal):
        # Define inclination as 0 through 90
        test_incl = proposal['value']
        if ((test_incl < 0) or (test_incl > 90)):
            raise traitlets.TraitError('incl must be in range 0 to 90')
        return test_incl

    @traitlets.observe('mass1', 'mass2', 'rad1', 'rad2', 'a', 'e', 'phi', 'N')
    def _refresh_orbit_model(self, change):
        # Method to call when any variable affecting orbit is called.  Since
        # the orbit change would affect light curves and radial velocity
        # curves, those are checked as well.
        # Update orbit if set to update continuously, update model(s), then
        # indcate updated, if not, set attribute to indicate not up to date
        #
        # Radius is included here because it can affect if the stars collide.
        if ((self.continuous_update) and (change['new'] != change['old'])):
            self.set_orbital_info()
            # Update radial velocity curve (since orbit change affects it)
            if (self.radvel_info is not None):
                self.set_radvel_info()

            # Update light curve data (since orbit change affects it)
            if (self.lc_info is not None):
                self.set_lc_info()

            self.up_to_date = True
        else:
            self.up_to_date = False

    @traitlets.observe('rv_sys', 'incl')
    def _refresh_rv_model(self, change):
        # Method to call when variables affecting the radial velocity
        # model (and thus the light curve model) get called.
        if ((self.continuous_update) and (change['new'] != change['old'])):
            # Revise values of radvel_info and lc_info if set
            if (self.radvel_info is not None):
                self.set_radvel_info()
            if ((change['name'] != 'rv_sys') and (self.lc_info is not None)):
                self.set_lc_info()
            self.up_to_date = True
        else:
            self.up_to_date = False

    @traitlets.observe('temp1', 'temp2', 'Na', 'Ntheta')
    def _refresh_only_lc_models(self, change):
        # Method to call when any attribute affecting only the light curve
        # model is changed.
        if ((self.continuous_update) and (change['new'] != change['old'])):
            # Revise values of radvel_info and lc_info if set
            if (self.lc_info is not None):
                self.set_lc_info()
            self.up_to_date = True
        else:
            self.up_to_date = False

        # Increment the number of models run since temp affect
        # appearance of model
        self.mdl_counter += 1

    ##
    # Define Methods
    ##

    def force_update(self):
        """
        Force update (useful if you turn off continuous updates)
        """
        # Update orbital information
        self.set_orbital_info()

        # Update radial velocity curve (since orbit change affects it)
        if (self.radvel_info is not None):
            self.set_radvel_info()

        # Update light curve data (since orbit change affects it)
        if (self.lc_info is not None):
            self.set_lc_info()

        self.up_to_date = True

    def set_orbital_info(self):
        """
        Assigns values for various parameters based on orbital model.  Allows
        computational Method _compute_orbital_info to be exposed for use
        outside class.
        """

        (self.P,
         self.ap,
         self.aa,
         self.maxrad,
         self.orbit_info) = self._compute_orbital_info(self.mass1,
                                                       self.mass2,
                                                       self.a_AU,
                                                       self.e,
                                                       self.phi,
                                                       self.N)

        # Check for collision and set that flag
        if ((self.rad1+self.rad2)*R_Sun > self.ap*AU):
            self.collision = True
        else:
            self.collision = False

        # Increment the number of models run
        self.mdl_counter += 1

    def _compute_orbital_info(self, mass1, mass2, a, e, phi=0, N=1000):
        """
        Using the masses of the two stars, their semi-major axis, and the
        eccentricity of the orbit, compute and return the orbital period,
        perihelion, aphelion separation, maximum distance from center of mass,
        and the orbital paths of the two stars in (x,y) coordinates in the
        plane of the orbit.

        If the inclination of orbital plane to the plane of the sky is zero,
        the z axis points toward us and we are viewing the system 'top down'.
        As the inclination angle increases, the angle between the z axis and
        line of sight equals the inclination angle.  I assume the y axis
        remains in the plane of the sky, so only the x and z axes undergo
        projection to the line of sight.

        No collision detection is done here, so two stars could collide if
        the sum of their radii is less than the periastron distance.

        This code is inspired in large by by Carrol and Ostlie's TwoStars code,
        although significant structural changes were made to take advantage of
        numpy's array arithmetic.

        Parameters
        ----------
        mass1 : float
                 mass of star 1 in solar masses
        mass2 : float
                 mass of star 2 in solar masses
        a : float
             semi-major axis of binary system reduced mass in AU
        e : float
             eccentricity of binary system reduced mass
        phi : float
               angle between semimajor axis and x axis in degrees
        N : int
             number of time steps to use in single period.

        Returns
        -------
        P : float
             orbital period in days
        ap : float
              Periastron distance in AU
        aa : float
              Apastron distance in AU
        maxrad : float
              maximum distance of either star from center of mass in AU
        OrbitInfo : Pandas dataframe
                     time series positions (in AU) and velocities (in km/s) of
                     both stars.
        """

        # Set up empty Pandas dataframe
        orbit_info = pd.DataFrame(columns=['time', 'r',
                                           'x1', 'y1',
                                           'vx1', 'vy1',
                                           'x2', 'y2',
                                           'vx2', 'vy2'])

        # Convert inputs into SI
        m1 = mass1*M_Sun
        m2 = mass2*M_Sun
        a_SI = a*AU
        phi_SI = phi*(np.pi/180.0)

        # Determine period using Kepler's Third Law (C&O 2.37)
        M = m1 + m2   # Total mass in system
        P = np.sqrt((4*np.pi*np.pi*a_SI*a_SI*a_SI)/(G*M))
        P_days = P/day   # Period in days

        # Compute the semimajor axes for each individual star's orbit
        mu = m1*m2/M     # Reduced mass (C&O 2.22)

        # Initiate the orbital computations
        L_ang = mu*np.sqrt(G*M*a_SI*(1 - e*e))  # Orbital ang. mom. (C&O 2.30)
        dAdt = L_ang/(2*mu)                     # Kepler's 2nd law (C&O 2.32)

        # Times covering entire orbit, including endpoints
        t = np.linspace(0, P, N+1, endpoint=True)
        dt = P/N  # The time step implemented above
        theta = np.zeros_like(t)
        r_SI = np.zeros_like(t)    # Empty positions array
        t_days = t/day  # Times in days

        # Loop through times to compute the radial position versus phase angle
        theta[0] = 0   # Initialize orbital phase angle
        for i in range(len(t)):
            r_SI[i] = a_SI*(1 - e*e)/(1 + e*np.cos(theta[i]))   # (C&O 2.3)
            if i < (len(t)-1):
                theta[i+1] = theta[i] + (2*dAdt/(r_SI[i]*r_SI[i]))*dt

        #
        # The rest of the time step parameters can be done in arrays.
        #

        # Velocity of reduced mass (C&O 2.3)
        v_SI = np.sqrt(G*M*(2/r_SI - 1/a_SI))
        v1_SI = (mu/m1)*v_SI
        v2_SI = (mu/m2)*v_SI

        # Position in orbital plane (x, y, z) frame [z=0]
        x_SI = r_SI*np.cos(theta + phi_SI)
        y_SI = r_SI*np.sin(theta + phi_SI)
        x1_SI = (mu/m1)*x_SI
        y1_SI = (mu/m1)*y_SI
        x2_SI = -(mu/m2)*x_SI
        y2_SI = -(mu/m2)*y_SI
        vx1_SI = -v1_SI*np.sin(theta + phi_SI)
        vy1_SI = v1_SI*np.cos(theta + phi_SI)
        vx2_SI = v2_SI*np.sin(theta + phi_SI)
        vy2_SI = -v2_SI*np.cos(theta + phi_SI)

        # Store in the data frame
        orbit_info['time'] = t_days
        orbit_info['r'] = r_SI/AU
        orbit_info['x1'] = x1_SI/AU
        orbit_info['y1'] = y1_SI/AU
        orbit_info['x2'] = x2_SI/AU
        orbit_info['y2'] = y2_SI/AU
        orbit_info['vx1'] = vx1_SI/1000
        orbit_info['vy1'] = vy1_SI/1000
        orbit_info['vx2'] = vx2_SI/1000
        orbit_info['vy2'] = vy2_SI/1000

        # Compute aphelion and perihelion
        ap_SI = np.min(r_SI)
        ap = ap_SI/AU
        aa_SI = np.max(r_SI)
        aa = aa_SI/AU

        # Compute farthest distance from origin by either object
        if (mass1 >= mass2):
            maxrad = (mu/m2)*aa
        else:
            maxrad = (mu/m1)*aa

        return (P_days, ap, aa, maxrad, orbit_info)

    def set_radvel_info(self):
        """
        Fills the radial velocity Pandas dataframe.  Allows
        computational Method _compute_orbital_info to be exposed for use
        outside class.
        """
        self.radvel_info = self._compute_radvel_info(self.orbit_info,
                                                     self.incl,
                                                     self.rv_sys)

    def wipe_radvel_info(self):
        """
        Wipes out the radial velocity Pandas dataframe, which turns off
        the computation of the data.
        """
        self.radvel_info = None

    def _compute_radvel_info(self, orbit_info, incl, rv_sys):
        """
        Computes the radial velocity curve of a binary star system with given
        parameters.  Returns the radial velocity curve as a Pandas dataframe
        and assuming positive radial velocities are radially AWAY from Earth.

        Parameters
        ----------
        OrbitInfo : Pandas dataframe
                     time series positions (in AU) and velocities (in km/s) of
                     both stars.
        incl : float
             inclination angle of z axis to line of sight (in degrees)
        rv_sys : float
             radial velocity of the center of mass of the system (in km/s)

        Sets
        -------
        self.radvel_info  : Pandas dataframe
                     time series radial velocities (in km/s) of both stars.
        """

        # Create empty dataframe
        RadVelInfo = pd.DataFrame(columns=['time', 'phase', 'v1r', 'v2r'])
        RadVelInfo['time'] = orbit_info['time']
        RadVelInfo['phase'] = orbit_info['time']/float(orbit_info['time'][-1:])

        # Compute radial velocities based on inclination angle
        proj = np.sin(incl*deg2rad)
        RadVelInfo['v1r'] = -proj * orbit_info['vx1'] + rv_sys
        RadVelInfo['v2r'] = -proj * orbit_info['vx2'] + rv_sys

        return RadVelInfo

    def set_lc_info(self):
        """
        Fills the light curve Pandas dataframe.  Allows
        computational Method _compute_orbital_info to be exposed for use
        outside class.
        """
        self.lc_info = self._compute_lc_info(self.orbit_info,
                                             self.incl,
                                             self.rad1,
                                             self.rad2,
                                             self.temp1,
                                             self.temp2,
                                             self.Na,
                                             self.Ntheta)

    def wipe_lc_info(self):
        """
        Wipes out the light curve Pandas dataframe, which turns off
        the computation of the data.
        """
        self.lc_info = None

    def _compute_lc_info(self, orbit_info, incl, rad1, rad2, temp1, temp2,
                         Na=100, Ntheta=360):
        """
        Computes the light curve of a binary star system with given parameters
        assuming some solar-like limb darkening.  Returns the light curve as
        a Pandas dataframe.

        WARNING: This function doesn't work properly if the two stars are
        going to collide during their orbit.

        Parameters
        ----------
        OrbitInfo : Pandas dataframe
                     time series positions (in AU) and velocities (in km/s) of
                     both stars.
        incl : float
             inclination angle of z axis to line of sight (in degrees)
        rad1 : float
             radius of star 1 (in solar radii)
        rad2 : float
             radius of star 2 (in solar radii)
        temp1 : float
             effective temperature of star 1 (in K)
        temp2 : float
             effective temperature of star 2 (in K)
        Na : int
              number of annuli to use in computing eclipse percentage
              (default of 100)
        Ntheta : int
                  number of angular steps to use in computing eclipse fraction
                  (default of 360)

        Returns
        -------
        lc_info : Pandas dataframe
                     time series percentage of total system flux visible.
        """

        # Create dataframe with 100% of flux visible at all times
        LightCurveInfo = pd.DataFrame(columns=['time', 'phase', 'F_norm',
                                               'in_front'])
        LightCurveInfo['time'] = orbit_info['time']
        LightCurveInfo['phase'] = (orbit_info['time'] /
                                   float(orbit_info['time'][-1:]))
        LightCurveInfo['F_norm'] = np.ones_like(orbit_info['time'])
        LightCurveInfo['in_front'] = np.zeros_like(orbit_info['time'])

        # Check if the stars are going to collide, if they are, return a zeroed
        # out light curve
        if (np.min(orbit_info['r']*AU/R_Sun) < (rad1 + rad2)):
            LightCurveInfo['F_norm'] = np.zeros_like(orbit_info['time'])
            return LightCurveInfo

        # Convert information into SI units
        rad1_SI = rad1*R_Sun
        rad2_SI = rad2*R_Sun
        x1_SI = orbit_info['x1']*AU
        y1_SI = orbit_info['y1']*AU
        x2_SI = orbit_info['x2']*AU
        y2_SI = orbit_info['y2']*AU

        # Transform positions into the 'sky' frame (x' toward observer, y'
        # parallel to y in orbital frame, z' perpendicular to both, so y' and
        # z' give projected position on sky).  Using the same coordinate setup
        # as Carroll and Ostlie's TwoStars code.
        cosi = np.cos(incl*deg2rad)
        sini = np.sin(incl*deg2rad)
        xp1 = x1_SI*sini
        yp1 = y1_SI
        zp1 = -x1_SI*cosi
        # xp2 = x2_SI*sini  # Never used
        yp2 = y2_SI
        zp2 = -x2_SI*cosi

        ##
        # Estimate uneclipsed flux and store flux in each annulus
        ##
        # Compute inner radii and radial stepsize for annuli
        r1, dr1 = np.linspace(0, rad1_SI, Na+1, endpoint=True, retstep=True)
        r2, dr2 = np.linspace(0, rad2_SI, Na+1, endpoint=True, retstep=True)
        # Convert inner radii to central radii and drop last radius
        r1 = r1[0:-1]+(dr1/2.)
        r2 = r2[0:-1]+(dr2/2.)
        # Compute area of each annulus
        dA1 = 2*np.pi*r1*dr1
        dA2 = 2*np.pi*r2*dr2
        # Compute flux contributed by each annulus (store for future use)
        dF1 = dA1*self._FluxVRad(temp1, r1, rad1_SI)
        dF2 = dA2*self._FluxVRad(temp2, r2, rad2_SI)
        # Total flux in system
        F1 = np.sum(dF1)
        F2 = np.sum(dF2)
        F_tot = F1 + F2

        # Compute center to center distances in sky frame
        d = np.sqrt((yp2-yp1)**2 + (zp2-zp1)**2)

        # Compute eclipsed flux at eclipsed times
        theta = np.linspace(0, 2*np.pi, Ntheta, endpoint=False)
        timesteps_to_check = np.where(d <= rad1_SI + rad2_SI)[0]
        for t_idx in timesteps_to_check:
            # Save separation of these two stars (in plane of sky)
            dist = d[t_idx]

            # Set up initial parameters of the computation
            if (xp1[t_idx] > 0):  # Star 1 in front, star 2 eclipsed
                F = F1  # Start with flux of star 1
                # Get information on flux distribution in star 2
                r = r2
                dF = dF2
                Rf, Rb = rad1_SI, rad2_SI   # Radii of stars in front and back
                # Hold positions of foreground and background stars
                ypf, zpf = yp1[t_idx], zp1[t_idx]
                ypb, zpb = yp2[t_idx], zp2[t_idx]
                in_front = 1
            else:  # Star 2 in front, star 1 eclipsed
                F = F2  # Start with flux of star 2
                # Get information on flux distribution in star 1
                r = r1
                dF = dF1
                Rf, Rb = rad2_SI, rad1_SI   # Radii of stars in front and back
                # Hold positions of foreground and background stars
                ypf, zpf = yp2[t_idx], zp2[t_idx]
                ypb, zpb = yp1[t_idx], zp1[t_idx]
                in_front = 2

            # If completely eclipsed, save flux of foreground star and skip
            # the computation
            if (dist + Rb < Rf):
                LightCurveInfo['F_norm'][t_idx] = F / F_tot
                LightCurveInfo['in_front'][t_idx] = in_front
                continue

            # Determine radii for computations of obscured background disk
            if (dist < Rb - Rf):  # Foregnd disk entirely within backgnd disk
                r_outer = dist + Rf
                r_inner = dist - Rf
                if (r_inner < 0):
                    r_inner = 0
            else:    # Foreground disk covers part of background disk
                r_outer = Rb
                r_inner = 0

            # Determine fraction of each annulus of each background star that
            # is behind foreground star
            visible = np.ones_like(r)
            # indicies to double check for this star
            check_idx = np.where((r >= r_inner) & (r <= r_outer))[0]

            for annulus in check_idx:
                # Compute visible fraction of this annulus by checking
                # how many points on this annulus are outside foreground star
                # radius.
                yp_annulus = r[annulus]*np.cos(theta) + ypb
                zp_annulus = r[annulus]*np.sin(theta) + zpb
                uneclipsed = np.sqrt((yp_annulus-ypf)**2 +
                                     (zp_annulus-zpf)**2) > Rf
                visible[annulus] = len(theta[uneclipsed]) / len(theta)

            # Add the visible flux from the background star and save the flux
            F += np.sum(dF*visible)
            LightCurveInfo['F_norm'][t_idx] = F / F_tot
            LightCurveInfo['in_front'][t_idx] = in_front

        return LightCurveInfo

    def _FluxVRad(self, T, r, Radius):
        """
        This function computes the flux at a given radial distance r from the
        center of a star with total radius Radius.

        It accounts for Limb Darkening using the 'solar=like' model 109 of
        Van Hamme (1993).  Probably not quite right for other stars, but better
        than no limb darkening correction at all.

        Reference:
           Van Hamme, W. 1993, AJ, 106, 1096.

        Parameters
        ----------
        T : float
             temperature of star (in K)
        r : float or numpy array of floats
             distance from center of star (in m)
        Radius : float
                  Total radius of photosphere of star (in m)

        Returns
        -------
        S : float or numpy array of floats
             Flux at this radial distance
        """

        # Bolometric coefficients from Van Hamme paper, model 109
        # (Temperature 5750, log g = 4.5)
        # This will generally underestimate the correction for cooler stars.
        x = 0.648
        y = 0.207

        # Compute lumb darkening correction using logarithmic model
        mu = np.cos((r*np.pi)/(Radius*2))
        ld_corr = 1 - x*(1 - mu) - y*mu*np.log10(mu)

        # Return flux assuming blackbody emission
        return sigma * T**4 * ld_corr


##
# Define previous functions at top level using BinaryStarModel methods
##
def OrbitalInfo(mass1, mass2, a, e, phi=0, N=1000):
    """
    Returns orbital information for a particular binary star model, this
    function is a placeholder to keep old code working.

    Parameters
    ----------
    mass1 : float
             mass of star 1 in solar masses
    mass2 : float
             mass of star 2 in solar masses
    a : float
         semi-major axis of binary system reduced mass in AU
    e : float
         eccentricity of binary system reduced mass
    phi : float
           angle between semimajor axis and x axis in degrees
    N : int
         number of time steps to use in single period.

    Returns
    -------
    P : float
         orbital period in days
    ap : float
          Periastron distance in AU
    aa : float
          Apastron distance in AU
    maxrad : float
          maximum distance of either star from center of mass in AU
    OrbitInfo : Pandas dataframe
                 time series positions (in AU) and velocities (in km/s) of
                 both stars.
    """
    bsm = BinaryStarModel(mass1, mass2, a, e, phi=0, N=1000, a_in_AU=True)
    return bsm.P, bsm.ap, bsm.aa, bsm.maxrad, bsm.orbit_info


def RadVelInfo(orbit_info, incl, rv_sys=0):
    """
    Returns a Pandas dataframe of radial velocities (assuming positive away
    from us) for a given binary star model.

    Parameters
    ----------
    OrbitInfo : Pandas dataframe
                 time series positions (in AU) and velocities (in km/s) of
                 both stars.
    incl : float
         inclination angle of z axis of system to line of sight (in degrees)
    rv_sys : float
         radial velocity of the center of mass of the system (in km/s)

    Returns
    -------
    RadVelInfo : Pandas dataframe
                 time series radial velocities (in km/s) of both stars.
    """

    # Create a binary star model just for the purpose of calling the radial
    # velocity fit with our orbital information.
    bsm = BinaryStarModel()
    return(bsm._compute_radvel_info(orbit_info, incl, rv_sys))


def LightCurveInfo(orbit_info, incl, rad1, rad2, temp1, temp2,
                   Na=100, Ntheta=360):
    """
    Returns a Pandas dataframe of light curve information for a given binary
    star model.

    WARNING: This function doesn't work properly if the two stars are going to
    collide during their orbit.

    Parameters
    ----------
    OrbitInfo : Pandas dataframe
                 time series positions (in AU) and velocities (in km/s) of
                 both stars.
    incl : float
         inclination angle of z axis of system to line of sight (in degrees)
    rad1 : float
         radius of star 1 (in solar radii)
    rad2 : float
         radius of star 2 (in solar radii)
    temp1 : float
         effective temperature of star 1 (in K)
    temp2 : float
         effective temperature of star 2 (in K)
    Na : int
          number of annuli to use in computing eclipse percentage
          (default of 100)
    Ntheta : int
              number of angular steps to use in computing eclipse percentage
              (default of 360)

    Returns
    -------
    LightCurveInfo : Pandas dataframe
                 time series percentage of total system flux visible.
    """

    # Create a binary star model just for the purpose of calling the light
    # curve fit with our orbital information.
    bsm = BinaryStarModel()
    return(bsm._compute_lc_info(orbit_info, incl, rad1, rad2,
                                temp1, temp2, Na, Ntheta))


class BinaryStarViewer(traitlets.HasTraits):
    """
    A class to display a binary star model using a pythreejs Renderer
    object.  You pass it the BinaryStarModel object to configure the model
    displayed.

    Parameters
    ----------
    bsm : BinaryStarModel Instance
                 instance of a BinaryStarModel
    t_idx0 : int
         Initial time index
    draw_grid : boolean
         Indicates whether orbital plane grid should be drawn or not
    draw_orbits : boolean
         Indicates whether orbits should be drawn or not
    view_width : int
         width of Render object in pixels (default: 600)
    view_height : int
         height of Render object in pixels (default: 600)
    lock_scale : Boolean
         If true, the scale is locked to render to scale and not increase size
         of objects.

    multiplier: float
         If this is greater than 1, the stars are being increased in size and
         not to scale.
    """

    # Define your traits here, as class attributes, using mdl_counter
    # to track when the binary star model says there is a change.
    t_idx = traitlets.Integer(allow_none=True)
    radius1 = traitlets.Float(allow_none=True)
    radius2 = traitlets.Float(allow_none=True)
    temp1 = traitlets.Float(allow_none=True)
    temp2 = traitlets.Float(allow_none=True)
    incl = traitlets.Float(allow_none=True)
    draw_grid = traitlets.Bool(allow_none=True)
    draw_orbits = traitlets.Bool(allow_none=True)
    lock_scale = traitlets.Bool(allow_none=True)
    mdl_counter = traitlets.Integer(allow_none=True)
    _initialized = False

    # Constants
    _AU2RSun = AU/R_Sun

    def __init__(self, bsm, t_idx0=0, draw_grid=True, draw_orbits=True,
                 lock_scale=False, enable_zoom=False,
                 view_width=600, view_height=600):
        ##
        # Need the orbital model be input as a widget.
        ##

        # initialize the class by calling the superclass
        super(BinaryStarViewer, self).__init__()

        # Initialize some variables
        self.t_idx = t_idx0
        self._view_factor = 1.25  # Place viewer this many times max distance
        self.minratio = 50    # Minimum ratio between size of objects and xmax
        self.draw_grid = draw_grid
        self.draw_orbits = draw_orbits
        self.lock_scale = lock_scale
        self.enable_zoom = enable_zoom
        self.multiplier = 1   # Increase size of the stars by this factor

        # Grab values from BinaryStarModel (bsm) and convert positions into
        # Solar radii and store those values in DataFrame
        self.bsm = bsm
        self._update_from_bsm()

        # Save initial radius to scale all other radii to this
        self._init_rad1 = bsm.rad1
        self._init_rad2 = bsm.rad2

        # Initialize a flat surface to contain orbital plane (accounting for
        # size of orbital plane + star radius)
        maxdist = self.bsm.maxrad*self._AU2RSun + max(self.radius1,
                                                      self.radius2)
        (self._xmax, self._grid_step) = self._grid_setup(maxdist,
                                                         update_pos=False)

        # Rescales the multiplier so that stars are rendered larger-than-scale
        # if view area becomes too large while avoiding having the stars
        # overlap
        if (self.lock_scale):
            self.multiplier = 1
        else:
            if ((self.minratio*min(self.radius1, self.radius2) < self._xmax)
                and (self.minratio*max(self.radius1, self.radius2) <
                     self.bsm.aa*self._AU2RSun)):
                self.multiplier = (self._xmax /
                                   (self.minratio *
                                    min(self.radius1, self.radius2)))
            else:
                self.multiplier = 1

        # Set up the scales based on initial radii
        sc1x = self.multiplier*self.radius1/self._init_rad1
        sc2x = self.multiplier*self.radius2/self._init_rad2
        self._scale1 = (sc1x, sc1x, sc1x)
        self._scale2 = (sc2x, sc2x, sc2x)

        # Create stars
        self._star1 = StarMesh(self.temp1, self.radius1, self._scale1,
                               [bsm.orbit_info['x1_RSun'][self.t_idx],
                                bsm.orbit_info['y1_RSun'][self.t_idx], 0])
        self._star2 = StarMesh(self.temp2, self.radius2, self._scale2,
                               [bsm.orbit_info['x2_RSun'][self.t_idx],
                                bsm.orbit_info['y2_RSun'][self.t_idx], 0])

        # Set up the viewing distance
        self._view_dist = self._view_factor*self._xmax

        # Generate the scene
        self._bsm_scene = self._create_scene()

        # Define initial viewing position
        rad_angle = self.incl*deg2rad
        self.init_position = (self._view_dist*np.sin(rad_angle),
                              0.0,
                              self._view_dist*np.cos(rad_angle))

        # Creates the camera so you can see stuff (with z-axis oriented up
        # as consistent with inclination to line of sight)
        # NOTE:
        # Using OrthographicCamera instead of PerspectiveCamera to show
        # view from Earth (far away) where both stars are really at the same
        # distance.
        self._starcam = p3j.OrthographicCamera(right=self._view_dist,
                                               left=-self._view_dist,
                                               top=self._view_dist,
                                               bottom=-self._view_dist,
                                               near=0, far=2*self._view_dist,
                                               position=self.init_position,
                                               up=[-1, 0, 0])

        # Makes a controller for the starcam camera looking toward the origin
        self._controller = p3j.OrbitControls(controlling=self._starcam,
                                             autoRotate=False,
                                             enableRotate=False,
                                             enableZoom=self.enable_zoom,
                                             target=[0, 0, 0])

        # creates the object that gets displayed to the screen
        self.renderer = p3j.Renderer(camera=self._starcam,
                                     scene=self._bsm_scene,
                                     controls=[self._controller],
                                     width=view_width,
                                     height=view_height)

        # Done initailizing object
        self._initialized = True

    def _update_from_bsm(self):
        '''
        This method will update values from the binary star model when
        called.
        '''
        self.temp1 = self.bsm.temp1
        self.temp2 = self.bsm.temp2
        self.radius1 = self.bsm.rad1
        self.radius2 = self.bsm.rad2
        self.incl = self.bsm.incl
        self.mdl_counter = self.bsm.mdl_counter  # No models computed yet
        self.orbit_info = self.bsm.orbit_info
        self.orbit_info['r_RSun'] = self.bsm.orbit_info['r']*self._AU2RSun
        self.orbit_info['x1_RSun'] = self.bsm.orbit_info['x1']*self._AU2RSun
        self.orbit_info['y1_RSun'] = self.bsm.orbit_info['y1']*self._AU2RSun
        self.orbit_info['x2_RSun'] = self.bsm.orbit_info['x2']*self._AU2RSun
        self.orbit_info['y2_RSun'] = self.bsm.orbit_info['y2']*self._AU2RSun

    def _grid_setup(self, maxdist, update_pos=True):
        '''
        Set the grid separation based on the aphelion distance.
        '''
        order = pow(10, np.floor(np.log10(maxdist)))
        coeff = maxdist/order

        # Round up the maximum extent and then choose about 5 grid step
        xmax = np.ceil(coeff)*order
        grid_step = xmax/5

        # Move the surface grid up to avoid rendering issues
        if ((update_pos) and (self.draw_grid)):
            self._surfgrid.position = (0.0, 0.0, xmax/300)

        return(xmax, grid_step)

    def _draw_orbits_RSun(self):
        """
        Take an orbit_info Pandas dataframe (modified to have positions in
        R_solar, and return two Line objects representing the two stars orbits.

        Assumes availability of orbital positions in solar radii units
        instead of AU.
        """

        # Determine number of rows of orbit_info, set color of each segment
        N = len(self.orbit_info)
        orbitcolor = ['yellow']*N

        # Define line material
        line_mat = p3j.LineBasicMaterial(linewidth=1,
                                         vertexColors='VertexColors')

        # Build Line object for star 1 orbit (to hover above plane)
        if (self.incl > 80):
            offset_factor = 1e-3
        else:
            offset_factor = 5e-3
        star1orbit = offset_factor*self._xmax*np.ones((N, 3))
        star1orbit[:, 0] = self.orbit_info['x1_RSun']
        star1orbit[:, 1] = self.orbit_info['y1_RSun']
        vertices1 = star1orbit.tolist()
        orbit1_geom = p3j.Geometry(vertices=vertices1, colors=orbitcolor)
        orbit1_line = p3j.Line(geometry=orbit1_geom, material=line_mat)

        # Build Line object for star 2 orbit (to hover above plane)
        star2orbit = offset_factor*self._xmax*np.ones((N, 3))
        star2orbit[:, 0] = self.orbit_info['x2_RSun']
        star2orbit[:, 1] = self.orbit_info['y2_RSun']
        vertices2 = star2orbit.tolist()
        orbit2_geom = p3j.Geometry(vertices=vertices2, colors=orbitcolor)
        orbit2_line = p3j.Line(geometry=orbit2_geom, material=line_mat)

        return (orbit1_line, orbit2_line)

    @traitlets.observe('mdl_counter')
    def _update_orbit(self, change):
        """
        This function is called if orbit_info of stars changes to update all
        the parameters accordingly.

        Only perform update if initialization has completed.
        """
        if (self._initialized):
            # Make class values match those in binary star model
            self._update_from_bsm()

            # Update surface grid and axes (could change if orbit changes)
            maxdist = self.bsm.maxrad*self._AU2RSun + max(self.radius1,
                                                          self.radius2)
            (self._xmax, self._grid_step) = self._grid_setup(maxdist)

            # Generate flat surface and grid for perspective
            if (self.draw_grid):
                self._surf_new, self._surfgrid_new = xyplane(self._xmax,
                                                             self._grid_step)
                self._surf.geometry = self._surf_new.geometry
                self._surfgrid.children = self._surfgrid_new.children

            # Grab values from BinaryStarModel (bsm) and convert positions into
            # Solar radii and store those values in DataFrame
            self.orbit_info['r_RSun'] = self.orbit_info['r']*self._AU2RSun
            self.orbit_info['x1_RSun'] = self.orbit_info['x1']*self._AU2RSun
            self.orbit_info['y1_RSun'] = self.orbit_info['y1']*self._AU2RSun
            self.orbit_info['x2_RSun'] = self.orbit_info['x2']*self._AU2RSun
            self.orbit_info['y2_RSun'] = self.orbit_info['y2']*self._AU2RSun

            # Update lines representing stars' orbits
            if (self.draw_orbits):
                (self._orbit1_new, self._orbit2_new) = self._draw_orbits_RSun()
                self._orbit1_line.geometry = self._orbit1_new.geometry
                self._orbit2_line.geometry = self._orbit2_new.geometry

            # Update position of stars
            self._update_time(change)

            # Update grid and star sizes and star colors if needed
            self._update_appearance(change)

    @traitlets.observe('t_idx')
    def _update_time(self, change):
        """
        This function handles just changes to the time index.

        Only perform update if entire object has been initialized.
        """
        if (self._initialized):
            self._star1.position = [self.orbit_info['x1_RSun'][self.t_idx],
                                    self.orbit_info['y1_RSun'][self.t_idx],
                                    0]
            self._star2.position = [self.orbit_info['x2_RSun'][self.t_idx],
                                    self.orbit_info['y2_RSun'][self.t_idx],
                                    0]

    @traitlets.observe('radius1', 'radius2', 'temp1', 'temp2', 'incl',
                       'lock_scale')
    def _update_appearance(self, change):
        # Update only if entire object has been initialized
        if (self._initialized):
            # Make class values match those in binary star model
            self._update_from_bsm()

            # Update surface grid and axes (could change if radius of stars
            # even if orbit doesn't change)
            maxdist = self.bsm.maxrad*self._AU2RSun + max(self.radius1,
                                                          self.radius2)
            (self._xmax, self._grid_step) = self._grid_setup(maxdist)

            # Update flat surface and grid for perspective (if needed)
            if (self.draw_grid):
                self._surf_new, self._surfgrid_new = xyplane(self._xmax,
                                                             self._grid_step)
                self._surf.geometry = self._surf_new.geometry
                self._surfgrid.children = self._surfgrid_new.children

            # Rescales the multiplier so that stars are rendered
            # larger-than-scale if view area becomes too large while avoiding
            # having the stars overlap
            if (self.lock_scale):
                self.multiplier = 1
            else:
                if ((self.minratio*min(self.radius1, self.radius2) < self._xmax)
                    and (self.minratio*max(self.radius1, self.radius2) <
                         self.bsm.aa*self._AU2RSun)):
                    self.multiplier = (self._xmax /
                                       (self.minratio *
                                        min(self.radius1, self.radius2)))
                else:
                    self.multiplier = 1

            # update scale of stars
            sc1x = self.multiplier*self.radius1/self._init_rad1
            sc2x = self.multiplier*self.radius2/self._init_rad2
            self._star1.scale = (sc1x, sc1x, sc1x)
            self._star2.scale = (sc2x, sc2x, sc2x)

            # Adjust star color
            hexcolor1 = tc.rgb2hex(tc.temp2rgb(self.temp1))[0]
            StarMeshColor(self._star1, hexcolor1)
            hexcolor2 = tc.rgb2hex(tc.temp2rgb(self.temp2))[0]
            StarMeshColor(self._star2, hexcolor2)

            # Move the camera (and adjust view if necessary)
            self._camera_update()

    def reset_grid_and_orbit(self):
        # This function is meant to be called from outside when the
        # flag on drawing the grid or orbit is changed.
        self.renderer.scene = self._create_scene()
 
    def _create_scene(self):
        # Generate flat surface and grid for perspective
        if (self.draw_grid):
            self._surf, self._surfgrid = xyplane(self._xmax, self._grid_step)

        # Construct lines representing stars' orbits
        if (self.draw_orbits):
            (self._orbit1_line, self._orbit2_line) = self._draw_orbits_RSun()

        # Makes the scene environment
        the_kids = [self._star1, self._star2]
        if (self.draw_grid):
            the_kids.extend([self._surfgrid, self._surf])
        if (self.draw_orbits):
            the_kids.extend([self._orbit1_line, self._orbit2_line])

        return p3j.Scene(children=the_kids, background='black')

    def _camera_update(self):
        '''
        Update the camera for adjustments in gridsize and position since last
        projection
        '''

        # Define initial viewing position (adjusted for inclination)
        self._view_dist = self._view_factor*self._xmax
        rad_angle = self.incl*deg2rad
        self._starcam.position = (self._view_dist*np.sin(rad_angle),
                                  0,
                                  self._view_dist*np.cos(rad_angle))

        # Define viewing region
        self._starcam.right = self._view_dist
        self._starcam.left = -self._view_dist
        self._starcam.top = self._view_dist
        self._starcam.bottom = -self._view_dist
        self._starcam.far = 2*self._view_dist

        # Tell OrbitControl to update since position changed
        self._controller.exec_three_obj_method('update',)
