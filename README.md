# Scallop

[![Build status](https://ci.appveyor.com/api/projects/status/j8uhgoibkeewd785?svg=true)](https://ci.appveyor.com/project/yakir12/scallop-jl) [![Build Status](https://travis-ci.org/yakir12/Scallop.jl.svg?branch=master)](https://travis-ci.org/yakir12/Scallop.jl)

[![Coverage Status](https://coveralls.io/repos/yakir12/Scallop.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/yakir12/Scallop.jl?branch=master) [![codecov.io](http://codecov.io/github/yakir12/Scallop.jl/coverage.svg?branch=master)](http://codecov.io/github/yakir12/Scallop.jl?branch=master)

Ray tracing for scallop eyes!
`Scallop.jl` allows for ray tracing scallop eyes using [`Rayden.jl`](https://github.com/yakir12/Rayden.jl).

## How to install
1. If you haven't already, install [Julia](https://julialang.org/downloads/) -> you should be able to launch it (some icon on the Desktop or some such)
2. Start Julia -> a Julia-terminal popped up
3. Copy: `Pkg.clone("git://github.com/yakir12/Rayden.jl.git"); Pkg.clone("git://github.com/yakir12/Scallop.jl.git")` and paste it in the newly opened Julia-terminal, press Enter
5. You can close the Julia-terminal after it's done running

## How to use
The main way you could use this package is by:

1) editing `variables.jl`.
2) running `main.jl`.
3) running a plotting function to see the results.

You can iterate through these three steps and see how the variables affect the optical function of the eye.

1) **`variables.jl`:** All the variables that control the outcome of the program are in the file [`variables.jl`](./src/variables.jl). This file (along with all the other sources files) is located in a folder on your computer, and you can find that folder by running the following command in the Julia-terminal window: `joinpath(Pkg.dir("Scallop"), "src")`. Open that file and edit the *quantities* of the variables in this file (e.g. from ` source_distance = 1m` to `source_distance = 12.3m`). While you can change the *units* of the variables as well (e.g. from ` source_distance = 1m` to `source_distance = 1cm`), make sure there is no space between the quantity and the unit (i.e. `1m` is correct and `1 m` is incorrect). Valid units are: `nm, μm, mm, cm, m, km, inch, ft, mi, rad, °`.
2) **`main.jl`:** Start a Julia-terminal, and run: `include(joinpath(Pkg.dir("Scallop"), "src", "main.jl"))`. 
3) **plot:** Once that is done running, you can copy-paste one or all of the following plotting functions:
    1) `ploteye()`: A schamtics of the scallop eye with some rays. This is useful to get a quick feel for the shape and function of the eye. 
    2) `plotpixels()`: Heat-maps of the PSF for each of the retinae. 
    3) `plotfwhm_aperture()`: A plot of the angular FWHM (in degrees) as a function of the pupil diameter (in microns). 
    4) `plotfwhm_distance()`: A plot of the angular FWHM (in degrees) as a function of the source distance (in cm). 
    5) `plotfwhm_aperture_distance()`: Heat-maps of the FWHM (in degrees) as a function of pupil diameter (in microns) and source distance (in cm).
 
## Optional arguments
All the plotting functions accept optional arguments. While you do not need to set these arguments (because they have default values), it might affect the accuracy and execution times greatly:

1) `ploteye(n_rays = x)`: Where `x` is the number of plotted rays.
2) `plotpixels(n_rays = n)`: Where `x` is the number of rays used to generate the PSF (more takes more time but results in a smoother heat-map).
3) `plotfwhm_aperture(n_data = x1, n_rays = x2, angle_step = x3)`: Where `x1` is the number of aperture-values used to generate the plot, `x2` is the number of rays used to calculate the FWHM, and `x3` is the angular resolution used to calculate the FWHM.
4) `plotfwhm_distance(n_data = x1, n_rays = x2, angle_step = x3)`: Where `x1` is the number of distances used to generate the plot, `x2` is the number of rays used to calculate the FWHM, and `x3` is the angular resolution used to calculate the FWHM.
5) `plotfwhm_aperture_distance(n_data = x1, n_rays = x2, angle_step = x3)`: Where `x1` is the number of aperture-values and distances used to generate the heat-maps, `x2` is the number of rays used to calculate the FWHM, and `x3` is the angular resolution used to calculate the FWHM.

Note that the angular resolution, `angle_step`, has to be denoted with its units (e.g. `1°` or `1rad`). All these values have reasonable default values, a balance between accuracy and execution times. 
