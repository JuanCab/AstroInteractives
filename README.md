Astro Interactives
==================

[*Astro Interactives*](https://juancab.github.io/AstroInteractives/) is a series of Jupyter Notebooks meant to provide 
interactive demos and simulations for introductory astronomy courses.  These notebooks were initially developed by 
- Juan Cabanela (Professor for Physics and Astronomy at Minnesota State University Moorhead)
- Andrew Louwagie Gordon (Student at Minnesota State University Moorhead)
- Sam Holen (Student at Minnesota State University Moorhead)

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


Dependencies
------------

These Jupyter notebooks require a `jupyter` server installation as well as the following python packages:

- `ipywidgets` (version >= 7.2.0)
- `bqplot`
- `pythreejs` (version >= 1.0.0)
- `appmode` (This allows the Jupyter notebook to run as an app)

which in turn depend on

- `traitlets` (version >= 4.3.0)
- `traittypes`
- `numpy`
- `pandas`


Known Issues
------------
- In some interactives, the animations are buffered by the web browser so that they lag behind the user inputs.  Efforts have been made to optimize the code to avoid these problems.

Help / Documentation
--------------------

- The documentation for how the Python notebooks do what they do are in the code comments for now.  We hope to eventually add copies of some of the lab exercises we do so people can see these on context.
- Send us an email at cabanela@mnstate.edu if you need any help/information.
