# Copyright (C) 2017-2019 by Daniel Shapero <shapero@uw.edu>
# Modified by Nick Holschuh <nholschuh@amherst.edu> 06/28/2023
#          --- This was made for two reasons
#          * to include handling of 1D data for flowine modeling
#          * to enable the use of xarray data objects (which are easier to interpret)
#
# This file is part of icepack.
#
# icepack is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The full text of the license can be found in the file LICENSE in the
# icepack source directory or at <http://www.gnu.org/licenses/>.

r"""Functions for interpolating gridded remote sensing data sets to finite
element spaces"""

import numpy as np
import ufl
import firedrake
import xarray
from scipy.interpolate import RegularGridInterpolator


def _sample(dataset, X, method):
    
    ######### Many common polar datasets, if stored as a geotiff, must be
    ######### loaded with rioxarray. These typically have extraneous
    ######### coordinate variables, which we drop here. We first check to see
    ######### if there is an empty dimension, we index away that dimension, and
    ######### then we drop the extraneous coordinate variables.
       
    dataset = dataset.squeeze(drop=True)
    coord_drop_list = ['lat','lon','band','mapping','spatial_ref']
    
    for i in coord_drop_list:    
        try: dataset = dataset.drop(i)
        except: pass

    coord_opts = sorted(list(dataset.coords))   
    ########### We have to manage non-ascending dimensions
    for i in coord_opts:
        if dataset[i].values[0] > dataset[i].values[-1]:
            dataset = dataset.isel({i:slice(None, None, -1)})
    
    ########### Here we manage the 1-dimensional case
    if X.ndim == 1:
        i = coord_opts[0]
        dim_res = np.median(np.diff(dataset[i].values))
        dimmin = max(X.min() - 2 * dim_res, dataset[i].min().values)
        dimmax = min(X.max() + 2 * dim_res, dataset[i].max().values)
        
        dataset = dataset.sel({i:slice(dimmin,dimmax)})
    
    ########### Here we manage the 2-dimensional case
    elif X.ndim > 1:
        for ind0,i in enumerate(coord_opts):
            dim_res = np.median(np.diff(dataset[i].values))
            dimmin = max(X[:, ind0].min() - 2 * dim_res, dataset[i].min().values)
            dimmax = min(X[:, ind0].max() + 2 * dim_res, dataset[i].max().values)
            
            dataset = dataset.sel({i:slice(dimmin,dimmax)})
            
            
    ######### Here we assemble the coordinate information
    coord_list = []
    for ind0,i in enumerate(coord_opts):
        coord_list.append(dataset[i].values)
    
    ########## Here we make sure the coordinate data are in the proper order
    ########## to correspond with the order of dimensions of the data
    new_coord_opts = []
    new_coord_list = []
    coord_sizes = dataset.values.shape
    for ind0,i in enumerate(coord_sizes):
        for ind1,j in enumerate(coord_list):
            if len(j) == i:
                new_coord_opts.append(coord_opts[ind1])
                new_coord_list.append(j)    
    
    
    ########### Finally, we transpose the object so that the object has x
    ########### coordinates first and y coordinates second
    interpolater_dataset = dataset.values.T
    new_coord_list.reverse()
    new_coord_opts.reverse()
    
    interpolator = RegularGridInterpolator(new_coord_list, interpolater_dataset, method=method,
                                           bounds_error=0, fill_value=0)
    
    return interpolator(X, method=method)


def interpolate_xr(f, Q_new, method="linear"):
    r"""Interpolate an expression or a gridded data set to a function space

    Parameters
    ----------
    f : xarray dataarray or tuple of xarray dataarray
        The gridded data set for scalar fields or the tuple of gridded data
        sets for each component
    Q : firedrake.FunctionSpace
        The function space where the result will live

    Returns
    -------
    firedrake.Function
        A finite element function defined on `Q` with the same nodal values
        as the data `f`
    """
    if isinstance(f, (ufl.core.expr.Expr, firedrake.Function)):
        return firedrake.interpolate(f, Q)

    mesh = Q_new.mesh()
    element = Q_new.ufl_element()
        
    ################################### When you are interpolating onto the object,
    ################################### the number of points required depends on 
    ################################### the degree of the polynomials involved.    
    number_of_elements = Q_new.dim()
    polynomial_degree = np.array(Q_new.ufl_element().degree())  
    
    if len(element.sub_elements()) == 0:
        print('For extruded meshes -- start by interpolating your data to a the base mesh')
        print('and then simply interpolate your values onto the extruded mesh using firedrake.interpolate(mesh3d,V)')
        
    # Cannot take sub-elements if function is 3D scalar, otherwise shape will mismatch vertical basis
    # This attempts to distinguish if multiple subelements due to dimension or vector function
    if issubclass(type(element), firedrake.VectorElement):
        element = element.sub_elements()[0]

    V_new = firedrake.VectorFunctionSpace(mesh, element)
    try:
        X = firedrake.interpolate(mesh.coordinates, V_new).dat.data_ro[:,:2]
    except:
        X = firedrake.interpolate(mesh.coordinates, V_new).dat.data_ro[:]

    q = firedrake.Function(Q_new)

    if isinstance(f, xarray.DataArray):
        q.dat.data[:] = _sample(f, X, method)
    elif isinstance(f, tuple) and all(
        isinstance(fi, xarray.DataArray) for fi in f
    ):
        for i, fi in enumerate(f):
            q.dat.data[:, i] = _sample(fi, X, method)
    else:
        raise ValueError(
            "Argument must be an Xarray Dataarray or a tuple of dataarrays!"
        )

    #return q, X
    return q
