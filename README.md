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

### Try it online with [Binder](http://mybinder.org/)

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=index.ipynb)

### Available Notebooks

The currently available Jupyter Notebooks (in no particular order) are the following:

1. [Radioactive Decay Interactive](RadioactiveDecay/radioactive_decay.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FRadioactiveDecay%2Fradioactive_decay.ipynb))  - This activity is used to illustrate the concept of radioactive decay and half-life.
2. [Radioactive Isochrones Interactive](RadioactiveIsochrones/Isochrones.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FRadioactiveIsochrones%2FIsochrones.ipynb))  - This activity is used to illustrate how radioactive dataing works using radioactive isochrones.  It is meant to be used after introducing the students to the concept of half-lives.
3. [Center of Mass Interactive](Center_of_Mass/Center_of_Mass.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FCenter_of_Mass%2FCenter_of_Mass.ipynb)) - This simulation simply allows a student to adjust the mass of two stars and dynamically shifts them back and forth to illustrate the center of mass position.  **Note:** The star colors and relative radii are approximately accurate assuming they are main sequence stars!
4. [Luminosity Calculator Interactive](LuminosityCalculator/LuminosityCalculator.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FLuminosityCalculator%2FLuminosityCalculator.ipynb))  - This activity is used to illustrate to students how the radius and temperature of a star affect its luminosity.  **Note:** The colors of the stars accurately correspond to the temperature, but the stars produced might not have radii and temperatures that correspond to real stars.
5. [Blackbody Simulation](BlackbodySimulation/BlackbodySimulation.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FBlackbodySimulation%2FBlackbodySimulation.ipynb)) - This simulation shows a Planck spectrum with user controls for the temperature.  We eventually hope to add B-V color estimation, similar to the [Flash-based UNL Blackbody Curves activity](http://astro.unl.edu/classaction/animations/light/bbexplorer.html)
6. [HR Diagram Interactive](HRdiagram/Interactive_HR_Diagram.ipynb) ([Launch in AppMode](https://mybinder.org/v2/gh/JuanCab/AstroInteractives/master?urlpath=%2Fapps%2FHRdiagram%2FInteractive_HR_Diagram.ipynb)) - This set of interactive graphics allows a student to explore the data that goes into the HR diagram, then how main sequence fitting can be used to determine the distance to a star cluster, and finally, how detailed stellar evolution models can be used to also determine the age of a star cluster.  **Note:** This notebook uses all "real data" with the exception of the fact that to show a brightness diagram corresponding to an HR diagram, the distances (and corresponding brightnesses) were 'assigned.'


### Dependencies

These Jupyter notebooks require a `jupyter notebook` installation as well as the following python packages:

- `ipywidgets` (version >= 7.2.0)
- `bqplot`
- `pythreejs` (version >= 1.0.0)
- `appmode` (This allows the Jupyter notebook to run as an app)

which in turn depends on

- `traitlets` (version >= 4.3.0)
- `traittypes`
- `numpy`
- `pandas`


Known Issues
------------
- The tooltip in the radioactive decay interactive is not working properly.
- Display of equations can lag, allowing underlying LaTeX to be visible briefly.

Help / Documentation
--------------------

- The documentation for how the Python notebooks do what they do are in the code comments for now.  We hope to eventually add copies of some of the lab exercises we do so people can see these on context.
- Send us an email at cabanela@mnstate.edu if you need any help/information.
