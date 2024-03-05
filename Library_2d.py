from numpy import *
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from netCDF4 import Dataset
from pathlib import Path
import json
import scipy
import math
import time
import warnings
import os
import pandas as pd
import glob
import cartopy
from cartopy import config
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from skimage.measure import block_reduce
# Interpolation
from scipy.interpolate import interpn, RegularGridInterpolator
warnings.filterwarnings("ignore")

##########################################################################################################
# Conversion to Cartesian coordinates
def cartesian_toUTM(*args):
  #=======================================================================================================
  # This function has been adapted from the Matlab version implemented by François Beauducel
  # and available at the link https://it.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll
  #=======================================================================================================
  # Input arguments
  varargin = args
  nargin = len(varargin)
  if nargin < 1:
    raise Exception('Not enough input arguments.')
  if (nargin > 1): 
    lat = varargin[0]
    lon = varargin[1]
    v = 2
  if (np.size(lat)>1) and (np.size(lon)>1) and (np.size(lat) != np.size(lon)):
    raise Exception('LAT and LON must be the same size or scalars.')  
  #------------------------
  # Available datums
  #------------------------
  datums = {'wgs84': [6378137.0, 298.257223563],
      'nad83': [6378137.0, 298.257222101],
      'grs80': [6378137.0, 298.257222101],
      'nad27': [6378206.4, 294.978698214],
      'int24': [6378388.0, 297.000000000],
      'clk66': [6378206.4, 294.978698214]}
  df_models = pd.DataFrame(datums)
  #-------------------------
  # Constants for conversion
  #-------------------------
  # Conversion rad to deg
  D0 = 180/math.pi	
  # UTM scale factor
  K0 = 0.9996
  # UTM false East [m]
  X0 = 500000
  #-------------------------
  # By default we will refer to wgs84
  #-------------------------
  datum = 'wgs84'
  zone = []
  if (nargin==3):
    datum = str(varargin[2])
    ref = df_models.loc[:,df_models.columns == datum].to_numpy()
    A1 = float(ref[0])
    F1 = float(ref[1])
  else:
    datum = 'wgs84'
    ref = df_models.loc[:,df_models.columns == datum].to_numpy()
    A1 = float(ref[0])
    F1 = float(ref[1])
  # Phi = latitude (rad)
  p1 = lat/D0		
  # Lambda = longitude (rad)
  l1 = lon/D0		
  # UTM zone automatic setting
  F0 = np.round((l1*D0 + 183)/6)
  # Conversion
  B1 = A1*(1 - 1/F1)
  E1 = np.sqrt((A1*A1 - B1*B1)/(A1*A1))
  P0 = 0/D0
  # UTM origin longitude (rad)
  L0 = (6*F0 - 183)/D0
  # UTM false northern [m]
  Y0 = 1e7*(p1 < 0)
  N = K0*A1
  C = coef(E1,0)
  B = C[0]*P0 + C[1]*math.sin(2*P0) + C[2]*math.sin(4*P0) + C[3]*math.sin(6*P0) + C[4]*math.sin(8*P0)
  YS = Y0 - N*B
  C = coef(E1,2)
  L = np.log(np.tan(np.pi/4 + p1/2)*(((1 - E1*np.sin(p1))/(1 + E1*np.sin(p1)))**(E1/2)))
  z = np.arctan(np.sinh(L)/np.cos(l1 - L0)) +1j*np.log(np.tan(math.pi/4 + np.arcsin(np.sin(l1 - L0)/np.cosh(L))/2))
  Z = N*C[0]*z + N*(C[1]*np.sin(2*z) + C[2]*np.sin(4*z) + C[3]*np.sin(6*z) + C[4]*np.sin(8*z))
  xs = Z.imag + X0
  ys = Z.real + YS
  return xs, ys

def coef(e,m):
  #COEF Projection coefficients
  #	COEF(E,M) returns a vector of 5 coefficients from:
  #		E = first ellipsoid excentricity
  #		M = 0 for transverse mercator
  #		M = 1 for transverse mercator reverse coefficients
  #		M = 2 for merdian arc
  if m==0:
      c0 = np.array([-175/16384, 0, -5/256, 0, -3/64, 0, -1/4, 0, 1,
              -105/4096, 0, -45/1024, 0, -3/32, 0, -3/8, 0, 0,
              525/16384, 0,  45/1024, 0, 15/256, 0, 0, 0, 0,
              -175/12288, 0, -35/3072, 0, 0, 0, 0, 0, 0,
              315/131072, 0,  0, 0, 0, 0, 0, 0, 0]).reshape(5,9)
      
  elif m==1:
      c0 = np.array([-175/16384, 0, -5/256, 0,-3/64, 0, -1/4, 0, 1,
                1/61440, 0, 7/2048, 0, 1/48, 0, 1/8, 0, 0,
              559/368640, 0, 3/1280, 0, 1/768, 0, 0, 0, 0,
              283/430080, 0, 17/30720, 0, 0, 0, 0, 0, 0,
          4397/41287680, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(5,9)
  elif m==2:
      c0 = np.array([-175/16384, 0, -5/256, 0, -3/64, 0, -1/4, 0, 1,
            -901/184320, 0, -9/1024, 0, -1/96, 0, 1/8, 0, 0,
            -311/737280, 0, 17/5120, 0, 13/768, 0, 0, 0, 0,
              899/430080, 0, 61/15360, 0, 0, 0, 0, 0, 0,
          49561/41287680, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(5,9)
    
  c = np.zeros((len(c0), 1))
  for i in range(len(c0)):
      c[i] = np.polyval(c0[i,:],e)
  return c
########################################################################################
# Reading ".xyz" files 
def read_XYZ(filename):
  data = pd.read_csv(filename, engine = 'python', sep = '\s+', header = None)
  x = data.loc[:,0].to_numpy()
  y = data.loc[:,1].to_numpy()
  z = data.loc[:,2].to_numpy()
  return x, y, z
  
#===================================================================
# Spatial domain
#===================================================================

def domain(left_point, right_point, step_grid):
    #---------------------------------------------------------------
    # This function will compute a spatial 
    # domain with a certain dx
    #---------------------------------------------------------------
    # INPUTS:
        # left_point =  Initial point
        # right_point =  Final point
        # step_grid = Spacing between each point 
    # OUTPUTS:
        # x = spatial domain
    #---------------------------------------------------------------
    axis_domain = arange(left_point, right_point + step_grid, step_grid)
    return axis_domain

#==================================================================
# 2d Integrand function
#==================================================================

def integrand_2d(m, n, K, xp, yp, a, b, H):
    #--------------------------------------------------------------
    # This function defines the 2d integrand function as in
    # equation (16) from paper Nosov and Kolesov, 2011.
    # The integrand is evaluated in a single point xp in the 
    # spatial domain.
    #--------------------------------------------------------------
    # INPUTS:
        # m = discrete wavenumber along x (point in the integration domain)
        # n = discrete wavenumber along y (point in the integration domain)
        # K = upper bound for the support of the integral (truncation)
        # xp = discrete grid point along x (point in the cell)
        # yp = discrete grid point along y (point in the cell)
        # a = length of the cell
        # b = width of the cell
        # H = depth (assumed to be constant within the cell)
    # OUTPUTS:
        # fun_point = integrand at point (xp, yp)
    #--------------------------------------------------------------
    k = sqrt(m**2 + n**2)
    fun_point = (cos(m*K*xp)*sin(m*K*a)*cos(n*K*yp)*sin(n*K*b))/(m*n*K*cosh(k*K*H))
    return fun_point

#===================================================================
# Legendre-Gauss adaptive quadrature evaluated at one point 
#===================================================================

def gauss_2d(La, Lb, Lc, Ld, K, xp, yp, a, b, H):
    #---------------------------------------------------------------
    # This function will evaluate the 2d composite formula for the
    # Gauss-Legendre quadrature with n points at one single point 
    # (xp,yp) in the spatial domain. The support of the 
    # integral is partioned automatically and within each partition
    # the result will be computed. The integrals computed 
    # individually within each partion are finally summed up.
    #---------------------------------------------------------------
    # INPUTS:
        # La, Lb = right and left limit for integral support (m)
        # Lc, Ld = right and left limit for integral support (n)
        # K = upper bound for integration
        # (xp, yp) = spatial point in the domain where evaluating the integrand
        # a = length of the cell
        # b = width of the cell
        # H = depth (constant along the cell)
    # OUTPUTS:
        # gauss_point = evaluation of the composite formula  at point xp
    #---------------------------------------------------------------
    # Change of variables (the support must be [-1,1]x[-1,1])
    half1 = (Lb-La)/2
    mid1 = (La+Lb)/2
    #half2 = (Ld-Lc)/2
    #mid2 = (Lc+Ld)/2

    half2 = Ld / 2.
    mid2 = Lc[1:] - half2

    # Four points
    x_i = [-0.577350269189626, -0.577350269189626, 0.577350269189626,0.577350269189626]
    y_i = [-0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626]
    w_i = [1, 1, 1, 1]
  
    gauss_point = 0
    for j in range(len(w_i)):
         gauss_point += half1*half2*(w_i[j] * integrand_2d(half1*x_i[j] + mid1, half2*y_i[j] + mid2, K, xp, yp, a, b, H))
    return gauss_point

#======================================================================
# Complete integral at one point in the spatial domain (Taylor + G-L)
#======================================================================
def integration2d_point(a, b, xp, yp, B0, H):
    #------------------------------------------------------------------
    # This function will evaluate the 1d numerical integration at one 
    # point xp in the spatial domain as a sum of two integrals:
    #    1. The first will integrate the Taylor expansion of the 
    #       integrand around 0
    #    2. The second will be given by the Gauss-Legendre composite
    #       formula 
    #------------------------------------------------------------------
    # INPUTS:
        # a = length of the cell
        # b = width of the cell
        # (xp, yp) = spatial point in the domain where evaluating the integrand
        # B0 = sea-floor deformation
        # H = depth (constant along the cell)
    # OUTPUTS:
        # elevation_point = evaluation of the composite formula  at point xp
    #--------------------------------------------------------------------
    # Upper bounds (integrals supports):
    #-----------------
    # Taylor series:
    eps = 1e-9
    # Gauss-Legendre:
    K = 5/H
    # Second integral
    #-----------------
    # Number of partitions
    npc = 2
    np_m = max([npc*int(K*max([abs(xp), a])/(2*pi)), 10])  
    np_n = max([npc*int(K*max([abs(yp), b])/(2*pi)), 10])
    # Coefficient
    coeff = 4.0*B0/pi**2
    # Integral support
    sup_m = linspace(eps/K, 1.0, np_m + 1)# %nrows
    sup_n = linspace(eps/K, 1.0, np_n + 1)# %ncols
    dsup_n = (1 - eps/K)/np_n
    # Adaptive Gauss-Legendre
    int1 = 0
    for i in range(len(sup_m)-1):
      #for j in range(len(sup_n)-1):
        gauss_point = gauss_2d(sup_m[i], sup_m[i+1], sup_n, dsup_n, K, xp, yp, a, b, H)
        int1 = int1 + gauss_point
    #----------------
    # First integral
    #-----------------
    int0 = a*b*(eps/K)**2
    #----------------------------------------
    # Summing first and second integral
    elevation_point = K*coeff*(int0 + sum(int1))
    return elevation_point

#===================================================================
# 2d Integration over the whole spatial domain
#===================================================================

def integration_2d(a, b, x, y,  B0, H):
    #--------------------------------------------------------------
    # This function will evaluate the 2d numerical integration at 
    # every point (xp, yp) in the spatial domain as defined in
    # the function integration_point
    #---------------------------------------------------------------
    # INPUTS:
        # a = length of the cell
        # b = width of the cell
        # x = discretization of the spatial domain along x
        # y = discretization of the spatial domain along y
        # B0 = bottom deformation within the cell (constant along the cell)
        # H = depth (constant along the cell)
    # OUTPUTS:
        # elevation = evaluation of the composite formula at every point
        # time_ex = execution time
    #---------------------------------------------------------------
    elevation = zeros((len(x), len(y)))
    start_time = time.time()
    for yp in range(len(y)):
        for xp in range(len(x)):
            elevation[xp, yp] = integration2d_point(a, b, x[xp], y[yp], B0, H)
    time_ex = start_time - time.time()
    return elevation, time_ex


#==================================================================================
# 2d Superposition
#==================================================================================

def superposition_2d(spatialX_domain, spatialY_domain,  dx, dy,  ncellX, ncellY, n_cell,  a, b, B0, H):
    #---------------------------------------------------------------
    # This function will evaluate the 2d superposition over the 
    # whole spatial domain. Within each cell an extended local 
    # domain is constructed, depending on the cell size and on the
    # value of the local depth, assumed to be constant within the 
    # cell. The local domain is symmetric and centered in zero and 
    # its coordinates are not referenced to the spatial domain.
    # Then, the result is shifted to be correctly placed in the 
    # spatial domain and summed to the others.
    #---------------------------------------------------------------
    # INPUTS:
        # spatialX_domain, spatialY_domain = entire domain (over which bathymetry is defined) - along both axes
        # dx = spacing between grid points of the spatial domain along x direction
        # dy = spacing betweem grid points of the spatial domain along y direction
        # ncellX = number of cells (along x)
        # ncellY = number of cells (along y)
        # n_cell = number of points within each cell (equal for the two directions)
        # a = length of the cell
        # b = width of the cell
        # B0 = bottom deformation within the cell (constant along the cell)
        # H = depth (constant along the cell)
    # OUTPUTS:
        # filtered_elev = evaluation of the composite formula over the domain
    #----------------------------------------------------------------
    # Initializing the array for storing the filtered elevation
    filtered_elev = np.zeros((len(spatialX_domain), len(spatialY_domain)))

    for j in range(ncellY):
        for i in range(ncellX):
            print(i,j)
            if B0[i,j]>= 0.1*B0.max() or B0[i,j]<=0.1*B0.min() and H[i,j]>1e3:#0.25*max(a,b):
                # Extended local domain for the cell (symmetric, centered in zero)
            
                xcell = domain(-4*H[i,j] - max(a,b), 4*H[i,j] + max(a,b), dx)
                ycell = domain(-4*H[i,j] - max(a,b), 4*H[i,j] + max(a,b), dy)
                # Shifting the extended local domain to be not centered in zero
              
                x_shift = spatialX_domain[n_cell*i] 
                y_shift = spatialY_domain[n_cell*j] 
                xcell_shifted = xcell + x_shift
                ycell_shifted = ycell + y_shift
                # Finding the correct position in the total domain
              
                new_posX = n_cell*(i) - int(4*H[i,j]/dx)
                new_posY = n_cell*(j) - int(4*H[i,j]/dy)
                # Interpolation
              
                iX = spatialX_domain[new_posX - n_cell:new_posX + int((8*H[i,j] + 2*max(a,b))/dx) + n_cell]
                iY = spatialY_domain[new_posY - n_cell:new_posY + int((8*H[i,j] + 2*max(a,b))/dy) + n_cell]
              
                if len(iX)>1 and len(iY)>1:
                    print(xcell.shape)
                    # Integral solution in the local extended domain
                    elev_cell, _ = integration_2d(a, b, xcell, ycell, B0[i, j], H[i, j])
                    #f = interpolate.interp2d(xcell_shifted, ycell_shifted, elev_cell, kind='linear')
                    #interp_values = f(iX, iY).T


                    interp_values = RegularGridInterpolator((xcell_shifted, ycell_shifted),elev_cell,
                                                        bounds_error=False, fill_value=None)


                    XX, YY = np.meshgrid(iX, iY)
                    interp_values = interp_values((XX, YY)).T


                    # Final solution
                    filtered_elev[new_posX-n_cell:new_posX+int((8*H[i,j]+2*max(a,b))/dx)+n_cell, new_posY-n_cell:new_posY+int((8*H[i,j]+2*max(a,b))/dy)+n_cell] += interp_values


                else:
                    pass
    return filtered_elev
