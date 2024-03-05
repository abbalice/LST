from Library_2d import *

#-----------------------------------------------Inizialization of arguments--------------------------------------------
def local_parser():
    parser = argparse.ArgumentParser(description=sys.argv[0])
    parser.add_argument('--event_loc', default=None, required=True, help='Specify the location of the event')
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
    event_loc = local_opts.event_loc
    return str(event_loc)

#------------------------------------------------------------------------------------------------------------------------
event_loc = from_local_parser()

workdir = os.getcwd()
# Defining diretory for outputs of this code
outDir = os.path.join(workdir, event_loc, 'outputs_step1')
try:
  os.makedirs(outDir)
except:
  os.path.exists(outDir)

# Retrieve bathymetric data
pathBathy = os.path.join(workdir, event_loc, 'bathymetry')
file_bathy = glob.glob(os.path.join(pathBathy,"*.xyz"))[0]
# Retrieve bathymetry and (lon, lat) coordinates
lon, lat, h = read_XYZ(file_bathy)
# Each of the lon/lat values are repeated: get the unique values
# and reshape the bathy array accordingly
lon_unique = np.unique(lon) 
lat_unique = np.unique(lat) 

[X, Y] = np.meshgrid(lon_unique, lat_unique)
h_matrix = np.reshape(h, X.shape) 

# Change reference system for the bathymetry (sea-depth should be defined positive)
H = -h_matrix
# Here i close the lines for the test

# Conversion to cartesian coordinates
x, y = cartesian_toUTM(lat_unique, lon_unique)
ix = np.argsort(x)
iy = np.argsort(y)
x = x[ix]
y = y[iy]

# Defining cell sizes
a = np.diff(x).mean()
b = np.diff(y).mean()

a = a/2
b = b/2

 # Defining the number of points for discretizing the cell
n_cell = 5

# Spacing between grid points
dx = 2*a/n_cell
dy = 2*b/n_cell

# Spatial domain
spatialX_domain = np.arange(min(x), max(x)+dx, dx)
spatialY_domain = np.arange(min(y), max(y)+dy, dy)

# Number of cells for domain discretization
ncellX = int((max(spatialX_domain) - min(spatialX_domain))/(2*a))
ncellY = int((max(spatialY_domain) - min(spatialY_domain))/(2*b))

print('Number of cells along x {}, Number of cells along y {}'.format(ncellX, ncellY))

#--------
# Temporarily export the positive-defined sea-bed for step 2
df_H = pd.DataFrame(H)
df_H.to_csv(os.path.join(outDir, 'H.txt'), sep=' ', header=False, index=None)

import json
data = {'a': a, 'b': b, 'dx': dx, 'dy': dy, 'x': list(spatialX_domain),  'y':list(spatialY_domain), 'ncellX':ncellX, 'ncellY':ncellY}
json_string = json.dumps(data, sort_keys=True, indent=8)

with open(os.path.join(outDir,'set_up.json'), 'w') as outfile:
    outfile.write(json_string)


