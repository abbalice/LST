from Library_2d import *
import time
#-----------------------------------------------Inizialization of arguments--------------------------------------------

def local_parser():
    parser = argparse.ArgumentParser(description=sys.argv[0])
    parser.add_argument('--workdir', default=None, required=True, help='Specify the working directory')
    parser.add_argument('--event_loc', default=None, required=True, help='Specify the location of the event')
    parser.add_argument('--component', default=None, required=True, help='Specify which component of the coseismic deformation')
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
    component = local_opts.component
    return str(workdir), str(event_loc), str(component)

#------------------------------------------------------------------------------------------------------------------------
workdir, event_loc, component = from_local_parser()
#------------------------------------------------------------------------------------------------------------------------
# Specify paths for retrieving the cosesmic deformation and the database of unit cells

pathDef = glob.glob(os.path.join(workdir, event_loc, 'coseismic_deformation', component, "*.xyz"))[0]
pathCells = os.path.join(workdir, 'database2_{}'.format(event_loc))

# Retrieve rbd (Residual Bottom Deformation) and convert it to a matrix
lon, lat, rbd = read_XYZ(pathDef)
lon_unique = np.unique(lon) 
lat_unique = np.sort(np.unique(lat))[::-1] 

[X, Y] = np.meshgrid(lon_unique, lat_unique)
B0 = np.reshape(rbd, X.shape) 

plt.figure()
plt.contourf(lon_unique,lat_unique,B0, cmap='coolwarm')
plt.colorbar()
plt.show()

# Retrieve corrected bathymetric data and model set-up
inpDir = os.path.join(workdir, event_loc, 'outputs_step1/')
H = pd.read_csv(inpDir+'H.txt', engine='python', sep='\s+', header=None, index_col=False).to_numpy()

with open(inpDir+'set_up.json') as json_file:
    data = json.load(json_file)
a = data['a']
b = data['b']
dx = data['dx']
dy = data['dy']
spatialX_domain = np.array(data['x'])
spatialY_domain = np.array(data['y'])
xdim, ydim = len(spatialX_domain), len(spatialY_domain)
n_cell=5
#-----------------------------
# Inizialization of the final filtered sea-surface elevation
filtered_elev = np.zeros((xdim, ydim))

# iterate over files in output directory
files = Path(pathCells).glob('*')

start_time = time.time()

for fid in files:
    fid = str(fid)
    # Retrieve cell indices
    cell = fid.split('Cell_')[1]
    i  = int(cell.split('_')[0])
    j = int(cell.split('_')[1])
    if H[i,j]>0:#0.25*max(a,b): 
       fp=os.path.join(fid,'cell_{}_{}.nc'.format(i,j))
       nc = Dataset(fp)
       xcell = nc['x'][:]
       ycell = nc['y'][:]
       elev_cell = nc['filter'][:,:]
       # Shifting the extended local domain to be not centered in zero
       x_shift = spatialX_domain[n_cell*i]
       y_shift = spatialY_domain[n_cell*j]
       xcell_shifted = xcell + x_shift
       ycell_shifted = ycell + y_shift
       # Finding the correct position in the total domain
       new_posX = n_cell*i - round(4*H[i,j]/(2*dx)) 
       new_posY = n_cell*j - round(4*H[i,j]/(2*dy)) 
       # Interpolation
       iX = spatialX_domain[new_posX:new_posX + round((8*H[i,j] + 2*max(a,b))/dx)] 
       iY = spatialY_domain[new_posY:new_posY + round((8*H[i,j] + 2*max(a,b))/dy)]
       interp_values = RegularGridInterpolator((xcell_shifted, ycell_shifted),elev_cell,
                                               bounds_error=False, fill_value=None)

       XX, YY = np.meshgrid(iX, iY)
       interp_values = interp_values((XX, YY)).T
       # Superposition
       filtered_elev[new_posX:new_posX+round((8*H[i,j] + 2*max(a,b))/dx), new_posY:new_posY+round((8*H[i,j] + 2*max(a,b))/dy)] += B0[i,j]*(interp_values)

end_time = time.time()
print('Execution time {} = {}[min.]'.format(event_loc, np.round((end_time-start_time)/60, 2)))
# Defining the directory where to store the local database and plots
outDir = os.path.join(workdir, event_loc, 'initial_conditions')

try:
  os.makedirs(outDir)
except:
  os.path.exists(outDir)

xi = block_reduce(filtered_elev, block_size=(n_cell, n_cell), func=np.mean)
#from scipy.ndimage.filters import uniform_filter
#xi = uniform_filter(xi, (n_cell, n_cell))


h = -H
h[h < 0] = 0

# Plot

from matplotlib.colors import Normalize

fig = plt.figure(figsize=(11,9))
orig_map=plt.cm.get_cmap('coolwarm')

# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
# Coseismic Deformation
gs = fig.add_gridspec(1, 2)
img_extent = [lon_unique.min(), lon_unique.max(), lat_unique.min(), lat_unique.max()]
ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
ax1.set_extent(img_extent, crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', color='k')
bathy1 = ax1.contourf(lon_unique, lat_unique, h, 100, transform=cartopy.crs.PlateCarree(), cmap="gray")
if abs(B0.min()) > B0.max():
   plot1 = ax1.contourf(lon_unique, lat_unique, B0, transform=cartopy.crs.PlateCarree(), alpha = 0.5, cmap=reversed_map)
else:
   plot1 = ax1.pcolormesh(X, Y, B0, transform=cartopy.crs.PlateCarree(), alpha = 0.5, cmap=orig_map, vmin=B0.min(), vmax=B0.max(), shading='auto')
ax1.contour(lon_unique, lat_unique, B0, transform=cartopy.crs.PlateCarree(), alpha = 0.5)
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=.2, color="k", alpha=0.5, linestyle="--")
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.ylines = True
gl.xlines = True
ax1.set_title('Coseismic deformation \n Max = {}, Min = {}'.format(np.round(B0.max(),2), np.round(B0.min(),2)))
# Filtered free-surface
ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
ax2.set_extent(img_extent, crs=ccrs.PlateCarree())
ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', color='k')
bathy2 = ax2.contourf(lon_unique, lat_unique, h, 100, transform=cartopy.crs.PlateCarree(), cmap="gray")
if abs(B0.min()) > B0.max():
   plot2 = ax2.contourf(lon_unique, lat_unique, xi, transform=cartopy.crs.PlateCarree(), alpha = 0.5, cmap=reversed_map)
else:
   plot2 = ax2.pcolormesh(X, Y, xi, transform=cartopy.crs.PlateCarree(), alpha = 0.5, cmap=orig_map, vmin=B0.min(), vmax=B0.max(), shading='auto')
ax2.contour(lon_unique, lat_unique, xi, transform=cartopy.crs.PlateCarree(), alpha = 0.5)
gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=.2, color="k", alpha=0.5, linestyle="--")
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.ylines = True
gl.xlines = True
ax2.set_title('Filtered sea-surface \n Max = {}, Min = {}'.format(np.round(xi.max(),2), np.round(xi.min(),2)))

fig.subplots_adjust(wspace=0.2)
cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
cbar=fig.colorbar(plot1, cax=cbar_ax, orientation='horizontal', extend='both')

fig.savefig(os.path.join(outDir, 'superposition_{}_{}'.format(event_loc, component)))
fig.show()
# Storing the database
name_file = 'Xi0_'+ event_loc + '_' + component + '.nc'
ncfile = Dataset(os.path.join(outDir, name_file), mode='w', format='NETCDF4')
ncfile.title = 'Nosov filter'
ncfile.subtitle = "Smoothed free-surface elevation"

x_dim = ncfile.createDimension('lon', len(lon_unique))     # longitude axis
y_dim = ncfile.createDimension('lat', len(lat_unique))    # latitude axis


# Define two variables with the same names as dimensions,
# a conventional way to define "coordinate variables".
x = ncfile.createVariable('lon', np.float32, ('lon',))
x.units = 'degrees'
x.long_name = 'Longitude'
y = ncfile.createVariable('lat', np.float32, ('lat',))
y.units = 'degrees'
y.long_name = 'Latitude'
#
# Define a 3D variable to hold the data
temp = ncfile.createVariable('oic',np.float64,('lon','lat')) # note: unlimited dimension is leftmost
temp.units = 'meters' 

x[:] = lon_unique
y[:] = lat_unique
# Write the data.  This writes the whole 3D netCDF variable all at once.
temp[:,:] = xi  # Appends data along unlimited dimension
print(xi.min(), xi.max(), B0.max(), B0.min())
ncfile.close()
print("Dataset created")

