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
StarMesh(temp, rad=,  scale, pos):
    Returns a pythreejs Mesh object corresponding to a star.
xyplane(max_dist, grid_space):
    Returns a pythreejs Mesh and SurfaceGrid corresponding to the xy plane.
axes(max_dist):
    Returns a pythreejs Line object representing the XYZ axes.
    
Created on Sun Jun 17 10:43:24 2018

@author: Juan Cabanela
(although many routines were originally written by Sam Holen and
Andrew Louwagie Gordon and just consolidated here)
"""

import numpy as np
import pythreejs as p3j
import tempNcolor as tc


def Rad_calc(mass = 1):
    '''
    Determines radius of a star given the mass.
    Uses approximate fit from ZAMS line


    Parameters
    ----------
    mass : float
            Mass of star in solar masses (default 1.0)

    Returns
    -------
    radius : float
              radius of star in units of solar radius. 
    '''
    return 10 ** (0.0757*(np.log10(mass))**4 -
                  0.1348*(np.log10(mass))**3 -
                  0.1355*(np.log10(mass))**2 +
                  0.8546*np.log10(mass) - 0.0516)


def Temp_calc(mass = 1):
    '''
    Determines the approximate temperature of a
    star given its mass. Uses a fit aquired from Eker et al (2015).

    Parameters
    ----------
    mass : float 
            Mass of star in solar masses (default 1.0)

    Returns
    -------
    temperature : float
                   temperature of star in Kelvin
    '''
    return 10 ** (0.6287*np.log10(mass) + 3.7404)


def ConfigStar(mass = 1.0):
    '''
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
    '''

    # Determine approximate radius in solar radii for both stars.
    radius = Rad_calc(mass)

    # Determines the approximate temperature of each star.
    temp = Temp_calc(mass)

    # Use scalar temperature to estimate hexcolor appropriate to each star
    hexcolor = tc.rgb2hex(tc.temp2rgb(temp))[0]

    return (radius, temp, hexcolor)


def StarMesh(temp=5700, rad=1,  scale=(1, 1, 1), pos=[0, 0, 0]):
    '''
    This function creates a pythreejs object that represents a star using
    a texture based on public domain STEREO Heliographic map made with
    images taken of the Sun on Dec. 30, 2011.  Image downloaded from
       https://stereo.gsfc.nasa.gov/360blog/
    and had its brightness and resolution rescaled.

    Parameters
    ----------
    temp : float
            temperature of star in Kelvin (default 5700)
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
    '''

    # Check is position is a list
    if isinstance(pos, list):
        # Check if this is a list of 3 items
        if (len(pos) != 3):
            raise TypeError('pos passed to StarMesh must be list of 3 numbers')
        # Check that all the items in the list are numbers
        for temp in pos:
            try:
                i = float(temp)
            except ValueError:
                raise TypeError('ValueError: pos contains list item that is not a number.')
    else:
        raise TypeError('pos passed to StarMesh must be list of 3 numbers')
    
    # Check is scale is a tuple
    if isinstance(scale, tuple):
        if (len(scale) != 3):
            raise TypeError('scale passed to StarMesh must be tuple of 3 numbers')
    else:
        raise TypeError('scale must be a tuple')
    
    # Define color and texture of star surface
    hexcolor = tc.rgb2hex(tc.temp2rgb(temp))[0]
    StarTexture = p3j.ImageTexture(imageUri='images/sun_surface.jpg')

    # Create sphere using MeshBasicMaterial (which is unaffected by lighting)
    StarSurface = p3j.MeshBasicMaterial(color=hexcolor, map=StarTexture)

    StarGeom = p3j.SphereBufferGeometry(radius=rad, widthSegments=32,
                                        heightSegments=16)
    return p3j.Mesh(geometry=StarGeom, material=StarSurface,
                    position=pos, scale=scale)


def xyplane(max_dist, grid_space):
    '''
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
    '''

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
    surface_material = p3j.MeshBasicMaterial(color='darkslategrey',
                                             transparent=True, opacity=0.5)
    surf = p3j.Mesh(geometry=surf_g, material=surface_material)
    grid_material = p3j.LineBasicMaterial(color='white')

    # To avoid overlap, lift grid slightly above the plane
    surfgrid = p3j.SurfaceGrid(geometry=surf_g, material=grid_material,
                               position=[0, 0, 1e-2])

    return surf, surfgrid


def axes(max_dist):
    '''
    Generate X, Y, Z axes of length max_width in the form of a pythreejs 
    Line object.
    
    Parameters
    ----------
    max_dist : float
                maximum extent of grid from origin in each dimension 

    Returns
    -------
    axes : pythreejs.Line
            a pythreejs Line object representing the xyz axes.
    '''
    axes_geom = p3j.Geometry(vertices=[[0, 0, 0], [max_width, 0, 0],
                                       [0, 0, 0], [0, max_width, 0],
                                       [0, 0, 0], [0, 0, max_width]],
                             colors = ['white', 'white', 'white',
                                       'white', 'white', 'white'])

    return p3j.Line(geometry=axes_geom,
                    material=p3j.LineBasicMaterial(linewidth=1))
