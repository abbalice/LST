from Library_1d import *

#-----------------------------------------------Inizialization of arguments--------------------------------------------
def local_parser():
    parser = argparse.ArgumentParser(description=sys.argv[0])
    parser.add_argument('--event_loc', default=None, required=True, help='Specify the location of the event')
    parser.add_argument('--component', default=None, required=True, help='Component of the coseismic deformation')
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
    component = local_opts.component
    return str(event_loc), str(component)

#------------------------------------------------------------------------------------------------------------------------
event_loc, component = from_local_parser()

workdir = os.getcwd()
# Retrieve bathymetric data
pathDef = glob.glob(os.path.join(workdir, event_loc, 'coseismic_deformation', component, "*.xyz"))[0]
pathBathy = os.path.join(workdir, event_loc, 'bathymetry')
file_bathy = glob.glob(os.path.join(pathBathy,"*.xyz"))[0]

# Retrieve rbd (Residual Bottom Deformation) and convert it to a matrix
lon, lat, rbd = read_XYZ(pathDef)

# Retrieve bathymetry and (lon, lat) coordinates
_, _, h = read_XYZ(file_bathy)
# Each of the lon/lat values are repeated: get the unique values
# and reshape the bathy array accordingly

lon_unique = np.unique(lon)
lat_unique = np.sort(np.unique(lat))[::-1]

[X, Y] = np.meshgrid(lon_unique, lat_unique)
h_reshaped = np.reshape(h, X.shape)
rbd_reshaped = np.reshape(rbd, X.shape)

x_utm, y_utm = cartesian_toUTM(lat_unique, lon_unique)

# Change reference system for the bathymetry (sea-depth should be defined positive)
h_positive = -h_reshaped


# Grid step for the transect
grid_step = 0.1

xx, yy, uzO = X.flatten(), Y.flatten(), rbd_reshaped.flatten()
grid_uzO = pygmt.xyz2grd(x=xx, y=yy, z=uzO, spacing=(0.01666667, 0.01666667), region=[min(lon), max(lon), min(lat), max(lat)])

xx, yy, bathy = X.flatten(), Y.flatten(), h_positive.flatten()
grid_bathy = pygmt.xyz2grd(x=xx, y=yy, z=bathy, spacing=(0.01666667, 0.01666667), region=[min(lon), max(lon), min(lat), max(lat)])

if component == 'vertical_2006' or component=='modelA_2006' or component == 'modelB_2006':
  #lon0, lon1 = 153, 155
  #lat0, lat1 = 47.5, 46.5

  #lon0, lon1 = 154, 156
  #lat0, lat1 = 48, 47

  lon0, lon1 = 152.5, 154
  lat0, lat1 = 46.8, 46

elif component == 'vertical_2007NW' or component=='modelA_2007NW':
  #lon0, lat0 = 152.5, 45.5
  #lon1, lat1 = 154, 44.75

  #lon0, lat0 = 153, 47
  #lon1, lat1 = 155.5, 45.8

  lon0, lat0 = 155, 47.5
  lon1, lat1 = 156, 47
 #-------------------------------------------------
elif component == 'vertical_2007SE' or component=='modelA_2007SE':
  #lon0, lat0 = 154.5, 47
  #lon1, lat1 = 156, 46

  #lon0, lat0 = 153, 46
  #lon1, lat1 = 154.5, 45

  lon0, lat0 = 153, 47
  lon1, lat1 = 155.5, 45.5

                                      
x0, y0 = cartesian_toUTM(lat0,lon0)
x1, y1 = cartesian_toUTM(lat1, lon1)

track_B0 = transect_elevation(lon0, lat0, lon1, lat1, grid_step, grid_uzO)
track_bathy = transect_elevation(lon0, lat0, lon1, lat1, grid_step, grid_bathy)
x = track_B0.p*1e3


# Defining the number of points for discretizing the cell
n_cell = 10

# Cell size
a = 1852 #np.diff(x).mean()
a = a/2
# Spacing between grid points
dx = 2*a/n_cell

distance = np.hypot((x1-x0), (y1-y0))

# Spatial domain
spatial_domain = domain(0, distance, dx)

# Number of cells for domain discretization
ncell = round((max(spatial_domain) - min(spatial_domain))/(2*a))

print('Number of cells along x {}'.format(ncell))


track_B0 = interp(spatial_domain, x, track_B0.elevation)
track_H0 = interp(spatial_domain, x, track_bathy.elevation)

B0 = block_reduce(track_B0, block_size=n_cell, func=np.mean)
H0 = block_reduce(track_H0, block_size=n_cell, func=np.mean)

start_time = time.time()
filtered_elev = superposition(spatial_domain, dx, ncell, n_cell, a, B0, H0, event_loc, component, plotIter = False)
time_sup = time.time()-start_time
#---------------------

xi_1d = block_reduce(filtered_elev, block_size=n_cell, func=np.mean)
spatial_domain = block_reduce(spatial_domain, block_size=n_cell, func=np.mean)
spatial_domain[-1] = spatial_domain[-2] + np.diff(spatial_domain).mean()


fig = plt.figure(figsize=(12, 8))
spec = fig.add_gridspec(nrows=1, ncols=1, width_ratios=[1], height_ratios=[1])

ax1 = fig.add_subplot(spec[0,0])
rect = ax1.patch
rect.set_facecolor('whitesmoke')
plt.text(x = x0, y = y0, s = "A", fontsize="medium")
plt.text(x = x1, y = y1, s = "A", fontsize="medium")
ax1.plot(spatial_domain*1e-3, B0, color='black', linewidth=2, label='Okada')
ax1.plot(spatial_domain*1e-3,  xi_1d, color='tab:blue', linewidth=2, label='1d-LST')
ax1.grid('on')
ax1.set_xlabel(r'$x$ [km]', fontsize=20)
ax1.set_xlim(0, (spatial_domain*1e-3).max())
ax1.set_ylabel(r'$\xi_{0}(x)$ [m]', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
ax1.legend(bbox_to_anchor=(0, 1.05, 1, 0.2), loc="lower left",
              mode="expand", borderaxespad=0, ncol=3, fontsize=12)

fig.savefig(os.path.join(os.getcwd(), '1d_transect_{}_{}.png'.format(event_loc, component)))
plt.show()


