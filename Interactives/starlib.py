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
OrbitalInfo(mass1, mass2, a, e, phi, N):
    Returns orbital period, perihelion, aphelion, and Pandas dataframe time
    series information of orbits.
RadVelInfo(orbit_info, incl):
    Returns radial velocity Pandas dataframe time series information 
    corresponding to a given orbit_info Pandas dataframe and inclination angle.
LightCurveInfo(orbit_info, incl, rad1, rad2, temp1, temp2, Na, Ntheta)
    Returns a light curve Pandas dataframe time series information
    corresponding to given orbit_info, inclination, and stellar parameters.
    
Created on Sun Jun 17 10:43:24 2018

@author: Juan Cabanela
(although many routines were originally written by Sam Holen and
Andrew Louwagie Gordon and just consolidated here)
"""

import numpy as np
import pandas as pd
import pythreejs as p3j
import tempNcolor as tc

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
    return 10 ** ( 0.0757*(np.log10(mass))**4
                  -0.1348*(np.log10(mass))**3
                  -0.1355*(np.log10(mass))**2
                  +0.8546*np.log10(mass))


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
                  -0.1600*np.log10(mass)**3
                  +0.2300*np.log10(mass)**2
                  +0.6084*np.log10(mass) + 3.7617)


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


def StarMesh(temp=Te_Sun, rad=1,  scale=(1, 1, 1), pos=[0, 0, 0]):
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
    hexcolor = tc.rgb2hex(tc.temp2rgb(float(temp)))[0]
    StarTexture = p3j.ImageTexture(imageUri='images/sun_surface.jpg')

    # Create sphere using MeshBasicMaterial (which is unaffected by lighting)
    StarSurface = p3j.MeshBasicMaterial(color=hexcolor, map=StarTexture)

    StarGeom = p3j.SphereBufferGeometry(radius=rad, widthSegments=32,
                                        heightSegments=16)
    return p3j.Mesh(geometry=StarGeom, material=StarSurface,
                    position=pos, scale=scale)


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
    surface_material = p3j.MeshBasicMaterial(color='darkslategrey',
                                             transparent=True, opacity=0.5)
    surf = p3j.Mesh(geometry=surf_g, material=surface_material)
    grid_material = p3j.LineBasicMaterial(color='white')

    # To avoid overlap, lift grid slightly above the plane scaling
    # by size of grid
    surfgrid = p3j.SurfaceGrid(geometry=surf_g, material=grid_material,
                               position=[0, 0, 1e-2*max_dist])

    return surf, surfgrid


def axes(max_dist):
    """
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
    """
    axes_geom = p3j.Geometry(vertices=[[0, 0, 0], [max_dist, 0, 0],
                                       [0, 0, 0], [0, max_dist, 0],
                                       [0, 0, 0], [0, 0, max_dist]],
                             colors=['white', 'white', 'white',
                                     'white', 'white', 'white'])

    return p3j.Line(geometry=axes_geom,
                    material=p3j.LineBasicMaterial(linewidth=1,
                                                   vertexColors='VertexColors'))


def OrbitalInfo(mass1, mass2, a, e, phi=0, N=1000):
    """
    Using the masses of the two stars, their semi-major axis, and the
    eccentricity of the orbit, compute and return the orbital period,
    perihelion, aphelion separation, maximum distance from center of mass,
    and the orbital paths of the two stars in (x,y) coordinates in the plane
    of the orbit.  
    
    If the inclination of orbital plane to the plane of the sky is zero, the 
    z axis points toward us and we are viewing the system 'top down'.  As the
    inclination angle increases, the angle between the z axis and line of sight
    equals the inclination angle.  I assume the y axis remains in the plane
    of the sky, so only the x and z axes undergo projection to the line
    of sight.
    
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
          Perihelion distance in AU
    aa : float
          Aphelion distance in AU
    maxrad : float
          maximum distance of either star from center of mass in AU
    OrbitInfo : Pandas dataframe 
                 time series positions (in AU) and velocities (in km/s) of 
                 both stars.                 
    """
    
    # Check the inputs
    if ((e<0) or (e>1)):
        raise ValueError('eccentricity of orbit must be between 0 and 1 (inclusive)')


    # Set up empty Pandas dataframe
    orbit_info = pd.DataFrame(columns=['time', 'r', 'x1', 'y1', 'vx1', 'vy1',
                                       'x2', 'y2', 'vx2', 'vy2'])

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
    L_ang = mu*np.sqrt(G*M*a_SI*(1 - e*e))   # Orbital angular mom. (C&O 2.30)
    dAdt = L_ang/(2*mu)                      # Kepler's 2nd law (C&O 2.32)

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
    orbit_info['r']  = r_SI/AU
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
    if (mass1>=mass2):
        maxrad = (mu/m2)*aa
    else:
        maxrad = (mu/m1)*aa

    return (P_days, ap, aa, maxrad, orbit_info)


def RadVelInfo(orbit_info, incl, rv_sys=0):    
    """
    Using the Pandas dataframe of orbital information (orbit_info) output
    by the OrbitalInfo().  It assumes the same coordinate system used for
    the orbital plane calculations in OrbitalInfo().  it assumes positive
    radial velocities are radially AWAY from Earth.

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
    
    # Create empty dataframe
    RadVelInfo = pd.DataFrame(columns=['time', 'phase', 'v1r', 'v2r'])
    RadVelInfo['time'] = orbit_info['time']
    RadVelInfo['phase'] = orbit_info['time']/float(orbit_info['time'][-1:])
   
    # Compute radial velocities based on inclination angle
    proj = np.sin(incl*deg2rad)
    RadVelInfo['v1r'] = -proj * orbit_info['vx1'] + rv_sys
    RadVelInfo['v2r'] = -proj * orbit_info['vx2'] + rv_sys
    
    return(RadVelInfo)



def LightCurveInfo(orbit_info, incl, rad1, rad2, temp1, temp2, Na=100, Ntheta=360):    
    """
    Using the Pandas dataframe of orbital information (orbit_info) output
    by the OrbitalInfo().  It assumes the same coordinate system used for
    the orbital plane calculations in OrbitalInfo().  it assumes positive
    radial velocities are radially AWAY from Earth.
    
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
    
    # Configure settings for computations
    Na = 100       # number of annuli to use to estimate flux
    Ntheta = 360   # Each annulus is broken into this many angular steps
        
    # Create dataframe with 100% of flux visible at all times
    LightCurveInfo = pd.DataFrame(columns=['time', 'phase', 'F_norm'])
    LightCurveInfo['time'] = orbit_info['time']
    LightCurveInfo['phase'] = orbit_info['time']/float(orbit_info['time'][-1:])
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
    # parallel to y in orbital frame, z' perpendicular to both, so y' and z' 
    # give projected position on sky).  Using the same frame as Carroll and 
    # Ostlie's TwoStars code.
    cosi = np.cos(incl*deg2rad)
    sini = np.sin(incl*deg2rad)
    xp1 = x1_SI*sini
    yp1 = y1_SI
    zp1 = -x1_SI*cosi
    xp2 = x2_SI*sini
    yp2 = y2_SI
    zp2 = -x2_SI*cosi

    ##
    ## Estimate uneclipsed flux and store flux in each annulus
    ##
    # Compute inner radii and radial stepsize for annuli
    r1, dr1 = np.linspace(0, rad1_SI, Na+1, endpoint=True, retstep = True)
    r2, dr2 = np.linspace(0, rad2_SI, Na+1, endpoint=True, retstep = True)
    # Convert inner radii to central radii and drop last radius
    r1 = r1[0:-1]+(dr1/2.)
    r2 = r2[0:-1]+(dr2/2.)
    # Compute area of each annulus
    dA1 = 2*np.pi*r1*dr1
    dA2 = 2*np.pi*r2*dr2
    # Compute flux contributed by each annulus (store for future use)
    dF1 = dA1*FluxVRad(temp1, r1, rad1_SI)
    dF2 = dA2*FluxVRad(temp2, r2, rad2_SI)
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
        if (xp1[t_idx]>0):  # Star 1 in front, star 2 eclipsed
            F = F1  # Start with flux of star 1
            # Get information on flux distribution in star 2
            r = r2
            dF = dF2
            Rf, Rb = rad1_SI, rad2_SI   # Radii of stars in front and back
            # Hold positions of foreground and background stars
            ypf, zpf = yp1[t_idx], zp1[t_idx]
            ypb, zpb = yp2[t_idx], zp2[t_idx]
        else:  # Star 2 in front, star 1 eclipsed
            F = F2  # Start with flux of star 2
            # Get information on flux distribution in star 1
            r = r1  
            dF = dF1
            Rf, Rb = rad2_SI, rad1_SI   # Radii of stars in front and back
            # Hold positions of foreground and background stars
            ypf, zpf = yp2[t_idx], zp2[t_idx]
            ypb, zpb = yp1[t_idx], zp1[t_idx]
            
        # If completely eclipsed, save flux of foreground star and skip 
        # the computation
        if (dist + Rb < Rf):
            LightCurveInfo['F_norm'][t_idx] = F / F_tot
            continue
        
        # Determine radii for computations of obscured background disk
        if (dist < Rb - Rf): # Foreground disk entirely within background disk
            r_outer = dist + Rf
            r_inner = dist - Rf
            if (r_inner < 0):
                r_inner = 0
        else:    # Foreground disk covers part of background disk
            r_outer = Rb
            r_inner = 0
        
        # Determine fraction of each annulus of each background star that is
        # behind foreground star
        visible = np.ones_like(r)
        # indicies to double check for this star
        check_idx = np.where((r>=r_inner) & (r<=r_outer))[0]
        
        for annulus in check_idx:
            # Compute visible fraction of this annulus by checking
            # how many points on this annulus are outside foreground star
            # radius.    
            yp_annulus = r[annulus]*np.cos(theta) + ypb
            zp_annulus = r[annulus]*np.sin(theta) + zpb
            uneclipsed = np.sqrt((yp_annulus-ypf)**2 + (zp_annulus-zpf)**2) > Rf
            visible[annulus] = len(theta[uneclipsed])/ len(theta)
            
        # Add the visible flux from the background star and save the flux
        F += np.sum(dF*visible)
        LightCurveInfo['F_norm'][t_idx] = F / F_tot

    return LightCurveInfo
    

def FluxVRad(T, r, Radius):
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
