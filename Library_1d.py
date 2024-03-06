from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import time
import warnings
import os
import pandas as pd
import utm
from skimage.measure import block_reduce
import pygmt
import glob
import sys
import argparse
from netCDF4 import Dataset
warnings.filterwarnings("ignore")

# Coordinates from transect (Nosov et al. 2008)
# max(lat) = 48.5, min(lon) = 152.5
# min(lat) = 46.2, max(lon) = 155.2
def cartesian_toUTM(*args):
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
  #     COEF(E,M) returns a vector of 5 coefficients from:
  #             E = first ellipsoid excentricity
  #             M = 0 for transverse mercator
  #             M = 1 for transverse mercator reverse coefficients
  #             M = 2 for merdian arc
  match m:
    case 0:
      c0 = np.array([-175/16384, 0, -5/256, 0, -3/64, 0, -1/4, 0, 1,
              -105/4096, 0, -45/1024, 0, -3/32, 0, -3/8, 0, 0,
              525/16384, 0,  45/1024, 0, 15/256, 0, 0, 0, 0,
              -175/12288, 0, -35/3072, 0, 0, 0, 0, 0, 0,
              315/131072, 0,  0, 0, 0, 0, 0, 0, 0]).reshape(5,9)

    case 1:
      c0 = np.array([-175/16384, 0, -5/256, 0,-3/64, 0, -1/4, 0, 1,
                1/61440, 0, 7/2048, 0, 1/48, 0, 1/8, 0, 0,
              559/368640, 0, 3/1280, 0, 1/768, 0, 0, 0, 0,
              283/430080, 0, 17/30720, 0, 0, 0, 0, 0, 0,
          4397/41287680, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(5,9)
    case 2:
      c0 = np.array([-175/16384, 0, -5/256, 0, -3/64, 0, -1/4, 0, 1,
            -901/184320, 0, -9/1024, 0, -1/96, 0, 1/8, 0, 0,
            -311/737280, 0, 17/5120, 0, 13/768, 0, 0, 0, 0,
              899/430080, 0, 61/15360, 0, 0, 0, 0, 0, 0,
          49561/41287680, 0, 0, 0, 0, 0, 0, 0, 0]).reshape(5,9)

  c = np.zeros((len(c0), 1))
  for i in range(len(c0)):
      c[i] = np.polyval(c0[i,:],e)
  return c

#----------------------------------------------------------------------------
# Read files with extension *.xyz
#----------------------------------------------------------------------------
def read_XYZ(filename):
  data = pd.read_csv(filename, engine = 'python', sep = '\s+', header = None)
  x = data.loc[:,0].to_numpy()
  y = data.loc[:,1].to_numpy()
  z = data.loc[:,2].to_numpy()
  return x, y, z

#------------------------------------------------------------------------------
# Get 1d Transect of a certain grid_file
#------------------------------------------------------------------------------
def transect_elevation(x0, y0, x1, y1, dx, grid_file):
  # Generate points along a great circle corresponding to the survey line
  # and store them in a pandas.DataFrame
  track_df = pygmt.project(
    unit = True,
    center="{}/{}".format(x0, y0),  # Start point of survey line (longitude/latitude)
    endpoint="{}/{}".format(x1, y1),  # End point of survey line (longitude/latitude)
    generate="{}".format(dx),  # Output data in steps of 0.1 degrees
    )
  track_df = pygmt.grdtrack(
    grid=grid_file,
    points=track_df,
    newcolname="elevation",
  )
  return track_df

#===================================================================
# 1d Spatial domain
#===================================================================

def domain(x_left, x_right, step_grid):
    #---------------------------------------------------------------
    # This function will compute a spatial 
    # domain with a certain dx
    #---------------------------------------------------------------
    # INPUTS:
        # x_left =  Initial point
        # x_right =  Final point
        # step_grid = Spacing between each point 
    # OUTPUTS:
        # x = spatial domain
    #---------------------------------------------------------------
    x = arange(x_left, x_right+step_grid, step_grid)
    return x

#==================================================================
# 1d Integrand function
#==================================================================

def integrand(m, K, xp, a, H):
    #--------------------------------------------------------------
    # This function defines the 1d integrand function as in
    # equation (8) from paper Nosov and Sementsov, 2014.
    # The integrand is evaluated in a single point xp in the 
    # spatial domain.
    #--------------------------------------------------------------
    # INPUTS:
        # m = discrete wavenumber (point in the integration domain)
        # K = upper bound for the support of the integral (truncation)
        # xp = discrete grid point (point in the cell)
        # a = length of the cell
        # H = depth (assumed to be constant within the cell)
    # OUTPUTS:
        # fun_xp = integrand at point xp
    #--------------------------------------------------------------
    fun_xp = (cos(m*K*xp)*sin(m*K*a))/(m*K*cosh(m*K*H))
    return fun_xp

#===================================================================
# Legendre-Gauss adaptive quadrature evaluated at one point 
#===================================================================

def gauss_1d(sup, dsup, K, n_gauss, xp, a, H):
    #---------------------------------------------------------------
    # This function will evaluate the 1d composite formula for the
    # Gauss-Legendre quadrature with n points at one single point 
    # xp in the spatial domain. The support of the 
    # integral is partioned automatically and within each partition
    # the result will be computed. The integrals computed 
    # individually within each partion are finally summed up.
    #---------------------------------------------------------------
    # INPUTS:
        # sup = discretized integral support
        # dsup = spacing between each point in sup
        # K = upper bound for the integral support (truncation)
        # n_gauss = number of points for the Gauss-Legendre quadrature
        # xp = spatial point in the domain where evaluating the integrand
        # a = length of the cell
        # H = depth (constant along the cell)
    # OUTPUTS:
        # gauss_xp = evaluation of the composite formula  at point xp
    #---------------------------------------------------------------
    # Change of variables (the support should be [-1,1])
    half = dsup/2.
    sup = sup[1:]-half 
    # Finding roots x_i of the legendre polynomials and weights w_i
    [x_i, w_i] = np.polynomial.legendre.leggauss(n_gauss)
    # Composite formula
    gauss_xp = 0
    for point in range(n_gauss):
        gauss_xp += half*(w_i[point] * integrand(sup + half*x_i[point], K, xp, a, H))
    return gauss_xp

#======================================================================
# Complete integral at one point in the spatial domain (Taylor + G-L)
#======================================================================
def integration_point(a, xp, B0, H):
    #------------------------------------------------------------------
    # This function will evaluate the 1d numerical integration at one 
    # point xp in the spatial domain as a sum of two integrals:
    #    1. The first will integrate the Taylor expansion of the 
    #       integrand around 0
    #    2. The second will be given by the Gauss-Legendre composite
    #       formula 
    #------------------------------------------------------------------
    # INPUTS:
        # sup = discretized integral support
        # dsup = spacing between each point in sup
        # xp = spatial point in the domain where evaluating the integrand
        # a = length of the cell
        # H = depth (constant along the cell)
    # OUTPUTS:
        # elevation_xp = evaluation of the composite formula  at point xp
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
    np = max([npc*int(K*max([abs(xp), a])/(2*pi)), 10])
    # Number of points for the Gauss-Legendre quadrature
    n_gauss = 3
    # Coefficient
    coeff = 2.0*B0/pi
    # Discretization of the integral support
    sup = linspace(eps/K, 1.0, np+1)
    # Spacing between points in the integral support
    dsup = (1 - eps/K)/np
    # Adaptive Gauss-Legendre
    int1 = gauss_1d(sup, dsup, K, n_gauss, xp, a, H)
    #----------------
    # First integral
    #-----------------
    int0 = a*eps/K
    #----------------------------------------
    # Summing first and second integral
    elevation_xp = K*coeff*(int0 + sum(int1))
    return elevation_xp

#===================================================================
# 1d Integration over the whole spatial domain
#===================================================================

def integration(a, x, B0, H):
    #--------------------------------------------------------------
    # This function will evaluate the 1d numerical integration at 
    # every point xp in the spatial domain as defined in
    # the function integration_point
    #---------------------------------------------------------------
    # INPUTS:
        # a = length of the cell
        # x = spatial domain
        # B0 = bottom deformation within the cell (constant along the cell)
        # H = depth (constant along the cell)
    # OUTPUTS:
        # elevation = evaluation of the composite formula at every point
        # time_ex = execution time
    #---------------------------------------------------------------
    elevation = zeros_like(x)
    start_time = time.time()
    for xp in range(len(x)):
        elevation[xp] = integration_point(a, x[xp], B0,  H)
    time_ex = start_time - time.time()
    return elevation, time_ex


#==================================================================================
# 1d Superposition
#==================================================================================

def superposition(spatial_domain, dx, ncell, nx_cell, a, B0, H, event_loc, component, plotIter = True):
    #---------------------------------------------------------------
    # This function will evaluate the 1d superposition over the 
    # whole spatial domain. Within each cell an extended local 
    # domain is constructed, depending on the cell size and on the
    # value of the local depth, assumed to be constant within the 
    # cell. The local domain is symmetric and centered in zero and 
    # its coordinates are not referenced to the spatial domain.
    # Then, the result is shifted to be correctly placed in the 
    # spatial domain and summed to the others.
    #---------------------------------------------------------------
    # INPUTS:
        # spatial_domain = entire domain (over which bathymetry is defined)
        # dx = spacing between grid points of the spatial domain
        # ncell = number of cells
        # nx_cell = number of points within each cell
        # a = length of the cell
        # B0 = bottom deformation within the cell (constant along the cell)
        # H = depth (constant along the cell)
        # plotIter = boolean to be set if one wants to plot at each iteration
        #            both the evaluation of the integral within the local 
        #            cell and its superposition over the spatial domain
    # OUTPUTS:
        # filtered_elev = evaluation of the composite formula over the domain
    #---------------------------------------------------------------

    # Initializing the array for storing the filtered elevation
    filtered_elev = zeros_like(spatial_domain)

    # Initializing the figure
    fig = plt.figure(figsize=(10,9))
    spec = fig.add_gridspec(nrows=2, ncols=1, width_ratios=[1], height_ratios=[1,1])
    ax = fig.add_subplot(spec[0,0])
    ax1 = fig.add_subplot(spec[1,0])

    start_time = time.time()
    for i in range(ncell):
        if H[i]>0:
            # Extended local domain for the cell (symmetric, centered in zero)
            xcell = domain(-4*H[i] - a, 4*H[i] + a, dx)
            # Solution to the integral within the cell
            elevn_cell, _ = integration(a, xcell, B0[i], H[i])
            x_domain = spatial_domain[nx_cell*i] + xcell
            # Finding the correct position in the total domain
            new_pos = nx_cell*i-round(4*H[i]/(2*dx))
            # Interpolation
            ix = spatial_domain[new_pos:new_pos + round((8*H[i] + 2*a)/dx)]
            interp_values = interp(ix, x_domain, elevn_cell)
            filtered_elev[new_pos:new_pos + round((8*H[i] + 2*a)/dx)] += interp_values
        #-------------------------------------------------------------------------------------------
        # Plotting at each iteration:
        if plotIter == True:
            #ax.plot(xcell, elevn_cell)    
            # The smoothed solution in the local extended domain
            ax.plot(x_domain, elevn_cell)
            ax.grid('on')
            #ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$\xi_{0}(x)$ [m]', fontsize=15)
            ax.set_title('Elevations within each cell (correctly placed)', fontsize=15)
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.tick_params(axis='both', which='minor', labelsize=15)
            # The result of the combination
            ax1.plot(spatial_domain, filtered_elev)
            ax1.grid('on')
            ax1.set_xlabel(r'$x$ [m]', fontsize=15)
            ax1.set_ylabel(r'$\xi_{0}(x)$ [m]', fontsize=15)
            ax1.set_title('Elevation smoothed', fontsize=15)
            ax1.tick_params(axis='both', which='major', labelsize=15)
            ax1.tick_params(axis='both', which='minor', labelsize=15)

    fig.savefig(os.path.join(os.getcwd(), 'Superposition_Iterations_{}_{}.png'.format(event_loc, component)))
    time_sup = time.time()-start_time
    #plt.show()
    return filtered_elev


def plot_grid(fig, spec, x0, y0, x1, y1, X, Y, U, bathy, component,event, profile_domain, profile_B0, profile_H):
    ax1 = fig.add_subplot(spec[0,0])
    plt.contourf(X, Y, U, 20, cmap = 'coolwarm', extend = 'both')
    #ax1.plot([x0, y0], [x1,y1], '--', color='black', linewidth=2)
    cbar = plt.colorbar();
    cbar.set_label('[m]')
    print(x0,y0,x1,y1)
    #ax1.plot([x0, y0], [x1,y1], 'o-', color='black', linewidth=2)
    ax1.set_xlabel('Lon [°], E')
    ax1.set_ylabel('Lat [°], N')
    ax1.set_title('Residual bottom deformation\n'+event)

    ax2 = fig.add_subplot(spec[0,1])
    plt.contourf(X, Y, bathy, 20, cmap='PRGn', extend = 'both')
    cbar = plt.colorbar();
    cbar.set_label('[m]')
    contours = plt.contour(X, Y, U, 10, colors='black')
    ax2.set_xlabel('Lon [°], E')
    ax2.set_ylabel('Lat [°], N')
    ax2.set_title('Bathymetry and contours\n'+event)

    ax3 = fig.add_subplot(spec[1,0])
    rect = ax3.patch
    rect.set_facecolor('whitesmoke')
    ax3.plot(profile_domain, profile_B0, color='black', linewidth=2)
    ax3.set_xlabel('x [km]')
    if component=='vertical':
        ax3.set_ylabel(r'$\eta_z$ [m]')
    ax3.set_title('Residual bottom deformation profile along transect')

    ax4 = fig.add_subplot(spec[1,1])
    rect = ax4.patch
    rect.set_facecolor('whitesmoke')
    ax4.plot(profile_domain, profile_H, color='black', linewidth=2)
    ax4.set_xlabel('x [km]')
    if component=='vertical':
        ax4.set_ylabel(r'$H$ [m]')
    ax4.set_title('Bathymetric profile along transect')

