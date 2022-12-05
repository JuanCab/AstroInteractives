Astro Interactives
==================

[*Astro Interactives*](https://juancab.github.io/AstroInteractives/) is a series of Jupyter Notebooks meant to provide 
interactive demos and simulations for introductory astronomy courses (and in a few cases, upper-division astrophysics courses).  These notebooks were initially developed by 
- Juan Cabanela (Professor for Physics and Astronomy at Minnesota State University Moorhead)
- Andrew Louwagie Gordon (Minnesota State University Moorhead Physics Major, Class of 2019)
- Sam Holen (Minnesota State University Moorhead Physics Major, Class of 2019)

This software is provide "as-is" and may contain bugs/errors that could  mis-represent astronomical reality.  You use this software at your own risk!

Goals
-----

-   Provide a freely downloadable set of Jupyter Notebooks that can be executed for use in introductory Astornomy labs and lecture activities to illustrate astronomical concepts.

Getting Started
---------------

### Try it online with [Binder](http://mybinder.org/) by clicking the button below

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?filepath=index.ipynb)


Available Notebooks
-------------------
The currently available Interactives are the following:

#### General Astronomy
1. [Small Angle Approximation Interactive](Interactives/SmallAngleEquation.ipynb) - This interactive was designed to allow students who are uncomfortable with mathematics visualize what is happening when the small angle approximation is valid.
2. [Flux vs.Luminosity Interactive](Interactives/FluxVsLuminositySimulation.ipynb) - This interactive is meant to be used in discussions of the concepts of flux/brightness, luminosity, and the inverse square law.
3. [Doppler Shift Interactives](Interactives/DopplerShift.ipynb) - These Doppler Shift interactives are designed to allow students to first model simple Doppler Shift and then to model the more complicated time-varying Doppler Shift seen in binary star and exoplanetary systems.

#### Solar System and Planetary Astronomy
1. [Radioactive Decay Interactives](Interactives/Radioactivity.ipynb)  - This pair of interactive plots is used to illustrate the concept of radioactive decay, half-life, and how these concepts can be used to determine the age of some rocks via the geochron method.
2. [Exoplanet System Simulation Interactive](Interactives/Exoplanet_Sim.ipynb) -  This interactive contains a full 3D simulation of a star system containing a single Jovian exoplanet, including a radial velocity or light curve for that system. 

#### Stellar Astronomy
1. [Center of Mass Interactive](Interactives/Center_of_Mass.ipynb) - This interactive model simply allows a student to adjust the mass and separation of two stars and see the shifts of their orbits, illustrating the center of mass position.
2. [Binary Star Simulation Interactive](Interactives/Binary_Star_Sim.ipynb) -  This interactive contains a full 3D simulation of a binary star system, including a radial velocity or light curve for that system.  The user can control most physical parameters of the system.
3. [Luminosity Calculator Interactive](Interactives/LuminosityCalculator.ipynb)  - This interactive is used to illustrate to students how the radius and temperature of a star affect its luminosity.  The colors of the stars accurately correspond to the temperature, but the stars produced might not have radii and temperatures that correspond to real stars.
4. [Blackbody Spectra Interactives](Interactives/Blackbody_Simulation.ipynb) - This set of interactive figures can be used to explore the properties of the Blackbody spectrum as applied to stars.  
5. [HR Diagram Interactives](Interactives/HR_Diagram.ipynb) - This set of interactive plots allows a student to explore the data that goes into the HR diagram, then how main sequence fitting can be used to determine the distance to a star cluster, and finally, how detailed stellar evolution models can be used to also determine the age of a star cluster.


#### Galactic and Extragalactic Astrophysics

The following interactives are meant to be used in an upper-division astrophysics class as a way of introducing l-v diagrams and 21-cm spectra.  These are not really meant to introduce the concepts but rather to allow students to explore how the distribution of neutral hydrogen gas both in position and velocity affects the observed l-v diagram of the Milky Way galaxy and the HI spectra of external galaxies.
1. [Synthetic l-v Diagram](Interactives/Synthetic_LV_Diagram.ipynb) -  This interactive takes a Milky Way-like rotation curve and neutral gas profile and generates a synthetic l-v diagram.  Users can then simply trace out a new rotation curve or neutral gas profile and see the corresponding l-v diagram.
2. [Synthetic HI Spectra](Interactives/Synthetic_Galaxy_HI_Spectra.ipynb) - This interactive model allows a student to see the single-dish (unresolved) HI spectra or the resolved HI spectra (aka velocity map) corresponding to a given model galaxy.  As with the l-v diagram interactive, users can  trace out a new rotation curve or neutral gas profile and see the corresponding spectra.


Dependencies
------------

These Jupyter notebooks require a `jupyter` server installation as well as the following python packages:

- `astropy`
- `bqplot`
- `pythreejs` (version >= 1.0.0)
- `appmode` (This allows the Jupyter notebook to run as an app within JupyterNotebook.  Another choice is to use `viola` which runs the Jupyter notebook as an app with its own server.)

which in turn depend on

- `ipywidgets` (version >= 8.0)
- `traitlets` (version >= 4.3.0)
- `traittypes`
- `numpy`
- `pandas`


Known Issues
------------
- In some interactives, the animations are buffered by the web browser so that they lag behind the user inputs.  Efforts have been made to optimize the code to avoid these problems.

Installation Instructions
-------------------------------
If you just want to run these notebooks locally on your own computer, you can use the following commands to do it:

1. Clone this repository with:
   `git clone https://github.com/JuanCab/AstroInteractives.git`
   or download the ZIP compressed copy of this repository if you don't have `git` and un-zip the repository into a new directory.

2. If you don't have a current python installation, install Anaconda Python 3 downloaded from [https://www.anaconda.com/download/](https://www.anaconda.com/download/)

3. Open a shell and add the `conda-forge` channel to the list of places conda looks for packages at using the command: 
   `conda config --append channel conda-forge`
    
4. Install the necessary packages not included with Anaconda (assuming you are using anaconda) with the command: 
   `conda install appmode astropy bqplot matplotlib pywidgets pythreejs` 
    or otherwise examine the `environment.yml` file to see the list of required python packages you will need to install.
    
5. Run the notebooks by changing to the directory containing this README file and typing: 
    `jupyter notebook index.ipynb`
    
6. **[Optional if you want to install these on a server for use by multiple users]** If you want to make these apps availabe to multiple students in a class, you may prefer to setup a Jupyterhub server.   We installed these files on a [The Littlest Jupyter Hub](https://github.com/jupyterhub/the-littlest-jupyterhub) installation on a virtual private server running Ubuntu linux.  I have [outlined one possible approach to doing this here](https://github.com/JuanCab/TLJH_AstroInteractives_Instructions).
Rather than create accounts for each students, we then use [Voila](https://github.com/QuantStack/voila) to serve the notebooks to our students without requiring them to login to the server.  It seems to work well.



Help / Documentation
--------------------

- The documentation for how the Python notebooks do what they do are in the code comments for now.  We hope to eventually add copies of some of the lab exercises we do so people can see these on context.
- Send us an email at cabanela@mnstate.edu if you need any help/information.
