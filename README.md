# petermann-icepack-flowmodel
Use icepack (Python finite element glacial flow model) to model dynamics of Petermann Glacier, NW Greenland



### Dependencies required for running the notebooks in this repository
icepack () and firedrake () -- install instructions are provided here: https://icepack.github.io/install/
Need up to date versions of: numpy, matplotlib, xarray, rioxarray, rasterio, geojson, sys, os, glob



### Additional run notes:
It is also important that you copy our interpolate_xr function into the icepack source directory, and modify init so it recognizes interpolate_xr as an available function.

