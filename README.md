# Computer simulations of incompressible fluids
Source code of my masters thesis as well as animated results. The goal of the thesis was to implement a Marker-And-Cell method algorithm based on Navier-Stokes equations. As the result program creates time simulations of incompressible fluids with various starting parameters.

## 📁 Code
Source code created in C++. Dividen into classes:
- ``System`` - main algorithm and settings,
- ``Type`` - enum type, describes types of cell on computational mesh,
- ``Cell`` - wrapper for ``Type``, cell in a mesh with flags,
- ``Particle`` - marked cell representing class,
- ``Matrix`` - additional class to represent matrices.

## 📁 Plots
Files that generate plots and animations in Gnuplot and Python (Matplotlib).

## 📁 Animations
Results for the configurations described in thesis in gif form. Some below as examples.

### *Lid-driven cavity* problem
![Lid-driven cavity problem animation](https://github.com/TMaczek/msc_thesis/blob/main/animations/6_2_anim.gif)

### Shallow wave breaking
![shallow wave breaking](https://github.com/TMaczek/msc_thesis/blob/main/animations/7_2_anim.gif)

### *Broken dam problem* with obstacle 
![broken dam with obstacle](https://github.com/TMaczek/msc_thesis/blob/main/animations/7_4_anim.gif)
