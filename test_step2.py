from Library_2d import *
from netCDF4 import Dataset
import pandas as pd
import json
import sys
import argparse

#-----------------------------------------------Inizialization of arguments--------------------------------------------

def local_parser():
    parser = argparse.ArgumentParser(description=sys.argv[0])
    parser.add_argument('--workdir', default=None, required=True, help='Working directory')
    parser.add_argument('--event_loc', default=None, required=True, help='Location of the database')
    parser.add_argument('--cellX', default=None, required=True, help='Cell in the domain along x')
    parser.add_argument('--cellY', default=None, required=True, help='Cell in the domain along y')
# load arguments
    args = parser.parse_args()
# if any
    if not sys.argv[1:]:
        print("Use -h or --help option for Help")
        sys.exit(0)
    return args
#-----------------------------------------
def from_local_parser():
    local_opts = local_parser()
    workdir = local_opts.workdir
    event_loc = local_opts.event_loc
    ncellX = local_opts.cellX
    ncellY = local_opts.cellY
    return str(workdir), str(event_loc), int(ncellX), int(ncellY)

workdir, event_loc, i, j = from_local_parser()
# Retrieve corrected bathymetric data and model set-up
inpDir = os.path.join(workdir, event_loc, 'outputs_step1/')
H = pd.read_csv(inpDir+'H.txt', engine='python', sep='\s+', header=None, index_col=False).to_numpy()
with open(inpDir+'set_up.json') as json_file:
    data = json.load(json_file)
a = data['a']
b = data['b']
dx = data['dx']
dy = data['dy']


# Defining the directory where to store the local database
outDir = 'database2_{}/Cell_{}_{}'.format(event_loc, i, j)
try:
  os.makedirs(outDir)
except:
  os.path.exists(outDir)

# The coseismic deformation is assumed scaled to one within each cell
B0 = 1

# Starting computation :  shallow depths will be avoided
if H[i,j]>0: #0.25*max(a,b):
  # Extended local domain for the cell (symmetric, centered in zero)
  xcell = np.arange(-4*H[i,j] - max(a,b), 4*H[i,j] + max(a,b) + dx, dx)
  ycell = np.arange(-4*H[i,j] - max(a,b), 4*H[i,j] + max(a,b) + dy, dy)
  # Evaluation of the integral within the cell c_ij
  elev_cell, _ = integration_2d(a, b, xcell, ycell, B0, H[i, j])
  # Store result within c_ij in a netcdf file
  ncfile = Dataset(os.path.join(outDir, 'cell_{}_{}.nc'.format(i, j)),mode='w',format='NETCDF4')
  ncfile.title='LSD'
  ncfile.subtitle="Sea-surface elevation in the LED"
  x_dim = ncfile.createDimension('x', len(xcell))     # longitude axis
  y_dim = ncfile.createDimension('y', len(ycell))    # latitude axis
  x = ncfile.createVariable('x', np.float32, ('x',))
  x.units = 'meters'
  x.long_name = 'Spatial domain along x axis'
  y = ncfile.createVariable('y', np.float32, ('y',))
  y.units = 'meters'
  y.long_name = 'Spatial domain along y axis'
  temp = ncfile.createVariable('filter',np.float64,('x','y'))
  temp.units = 'meters'
  x[:] = xcell
  y[:] = ycell
  temp[:,:] = elev_cell
  ncfile.close()
  print("Dataset created")
else:
  print("cell {} {} under threshold".format(i,j))
