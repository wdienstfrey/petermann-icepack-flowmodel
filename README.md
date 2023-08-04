# petermann-icepack-flowmodel
Uses icepack (Python finite element glacial flow model developed by Daniel Shapero) to model dynamics of Petermann Glacier, NW Greenland. Uses many of icepack's unique features such as the inverse solver, the 3D advection-diffusion heat transport forward model, and the hybrid flow model. Both the 2D flowline and the 3D basin model are defined by a monolayer extruded mesh. For questions about icepack refer to Daniel Shapero's github at https://github.com/icepack/icepack



### Dependencies required for running the notebooks in this repository
Need icepack () and firedrake () -- install instructions are provided here: https://icepack.github.io/install/

Up to date versions of: numpy, matplotlib, xarray, rioxarray, rasterio, geojson, sys, os, glob



### Data to run notebooks
All files except BedMachine can be found here: https://drive.google.com/drive/folders/1vSetbOF5Iy26ETMAF_-VXlrqotyp_3EO?usp=sharing

BedMachine and other file sources can be found here:

BedMachine: https://nsidc.org/data/idbmg4/versions/5

MEaSUREs: https://nsidc.org/grimp

Martos Geothermal heat flux: https://doi.pangaea.de/10.1594/PANGAEA.892973?format=html#download

RACMO2 Surface Temperature/Surface Mass Balance: https://doi.org/10.5194/tc-10-2361-2016 (Need to email author Brice Noel for dataset)

Hillshade: https://drive.google.com/drive/folders/1vSetbOF5Iy26ETMAF_-VXlrqotyp_3EO?usp=sharing (creator is Nicholas Holchuh)



### Additional run notes:
It is also important that you copy our interpolate_xr function into the icepack source directory, and modify init so it recognizes interpolate_xr as an available function.

