# petermann-icepack-flowmodel
Use icepack (Python finite element glacial flow model) to model dynamics of Petermann Glacier, NW Greenland



### Dependencies required for running the notebooks in this repository
icepack () and firedrake () -- install instructions are provided here: https://icepack.github.io/install/
\n numpy
\n matplotlib
\n xarray
\n rioxarray
\n rasterio
\n geojson
\n sys
\n os
\n glob



### Additional run notes:
It is also important that you copy our interpolate_xr function into the icepack source directory, and modify init so it recognizes interpolate_xr as an available function.

