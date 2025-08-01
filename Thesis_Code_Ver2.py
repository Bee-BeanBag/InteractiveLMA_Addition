#%% Imports
from __future__ import print_function
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.geodesic import Geodesic
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import csv
import datetime as dt
from datetime import datetime, timedelta
from geopy import distance
import h5py
from ipywidgets import interact, widgets
from math import radians, cos, sin, asin, sqrt, degrees, atan, log
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import matplotlib.dates as md
from matplotlib.dates import AutoDateLocator
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import Formatter, FormatStrFormatter, MaxNLocator
from metpy.plots import USCOUNTIES
import nexradaws
conn = nexradaws.NexradAwsInterface()
import numpy as np
import numpy.ma as ma
norm = np.linalg.norm
import os
import pandas as pd
from pathlib import Path
import pickle
from PIL import Image
import pyart
import pymap3d as pm
from pyxlma import coords
from pyxlma.lmalib.flash.cluster import cluster_flashes
import pyxlma.lmalib.flash.properties
from pyxlma.lmalib.flash.properties import flash_stats, filter_flashes
from pyxlma.lmalib.grid import  create_regular_grid, assign_regular_bins, events_to_grid
from pyxlma.lmalib.io import read as lma_read
from pyxlma.plot.interactive import InteractiveLMAPlot
from pyxlma.plot.xlma_base_plot import FractionalSecondFormatter
from pyxlma.plot.interactive import output as lma_plot_output
from pyxlma.plot.xlma_base_plot import subplot_labels, inset_view, BlankPlot
from pyxlma.plot.xlma_plot_feature import color_by_time, plot_points, setup_hist, plot_3d_grid, subset, plot_2d_network_points
import pyxlma.plot.xlma_super_plot_feature
import pyxlma.plot.xlma_super_base_plot
import pyxlma.plot.radar as lmarad
from pypalettes import load_cmap
import pyproj as proj4
import requests
import shapely.geometry as sgeom
from scipy import stats
from scipy import spatial
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import sys
import gc
from geodatasets import get_path
import geopandas
import glob
import gzip
import urllib.request
import warnings
warnings.filterwarnings("ignore")
import wradlib
import xarray as xr
import xradar as xd
import zipfile

try:
    from metpy.plots import USCOUNTIES
    county_scales = ['20m', '5m', '500k']
    COUNTIES = USCOUNTIES.with_scale(county_scales[0])
except ImportError:
    COUNTIES = None
        
#%% Run From Here

#Bee's
obsticles = "https://aeronav.faa.gov/Obst_Data/DOF_250511.zip"
urllib.request.urlretrieve(obsticles, "DOF_250511.zip")

zipdata = zipfile.ZipFile('DOF_250511.zip')
zipinfos = zipdata.infolist()

#Input ENTLN data here
entln_data = lma_read.entln('/Users/BenLa/.spyder-py3/LEEpulse_full.csv')

# iterate through each file
for zipinfo in zipinfos:
    # This will do the renaming
    if '.csv' in zipinfo.filename:
        print("Found ", zipinfo.filename)
    # We can override the filename. Neat!\n",
    zipinfo.filename = 'turbine_locations.csv'
    zipdata.extract(zipinfo)


def haversine(lat1, lon1, lat2, lon2):

      R = 6378.137e3 # this is in meters.  For Earth radius in kilometers use 6372.8 km

      dLat = radians(lat2 - lat1)
      dLon = radians(lon2 - lon1)
      lat1 = radians(lat1)
      lat2 = radians(lat2)

      a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
      c = 2*asin(sqrt(a))

      return R * c    

def negativeaxis(position, radar):
    
    if position < 0:
        if position < radar.longitude['data'][0]:
            return True
    else:
        if position < radar.latitude['data'][0]:
            return True
    return False

def explode_xy(xy):
    xl=[]
    yl=[]
    xy2=[]
    for y in range(len(xy)):
        xl.append(xy[y][0])
        yl.append(xy[y][1])
    xl = pm.geodetic2enu(yl, xl, 10000, yl[0], xl[0], 10000)[0]*1e-3
    yl = pm.geodetic2enu(yl, xl, 10000, yl[0], xl[0], 10000)[1]*1e-3
    for u in range(len(xl)): xy2.append([xl[u], yl[u]])
    return xy2

def ReadRadar(radar_file):
    radar_list = []
    start_times = []
    for files in radar_file:
        print("Reading: " + files)
        radar = pyart.io.read(files)
        radar_list.append(radar)
        radar_name = radar.metadata['instrument_name']
        start = dt.datetime.strptime(files.split('/')[-1], radar_name + "%Y%m%d_%H%M%S_V06").strftime('%d/%m/%Y %H%M.%S')
        start_times.append(start)
    
    return radar_list, start_times

def GridRadar(radar_objs):
        grided = []
        for radars in radar_objs:
                grided.append(pyart.map.grid_from_radars(radars, grid_shape=(20, 181, 181), grid_limits=((100, 15000), (-175000, 175000), (-155000, 155000)), grid_origin=(pyart.graph.RadarMapDisplay(radars).loc[0], pyart.graph.RadarMapDisplay(radars).loc[1]), fields=["reflectivity", 'spectrum_width', 'differential_reflectivity', 'cross_correlation_ratio', 'velocity', 'differential_phase']))
        
        return grided

#Bee's
data = zipfile.ZipFile('DOF_250511.zip')
infos = data.infolist()

for info in infos:
    if '36-NY' in info.filename:
        print(info.filename)
        info.filename = 'NewYork_Obsticles.Dat'
        data.extract(info)
        
for info in infos:
    if '46-SD' in info.filename:
        print(info.filename)
        info.filename = 'SouthDakota_Obsticles.Dat'
        data.extract(info)
        
#%%
#ds = xr.open_dataset('/Users/BenLa/.spyder-py3/FlashAnalysis/Flash_Files/LEE/LYLOUT_221117_080000_432000_map4000.nc')
#starttime, endtime = ds.grid_time_edge[0].data, ds.grid_time_edge[-1].data

#Imput LMA File Here 
ds = xr.open_dataset('/Users/BenLa/.spyder-py3/LYLOUT_221117_080000_360000_map4000NovFull1.nc')

starttime, endtime = ds.grid_time_edge[0].data, ds.grid_time_edge[-1].data

lma_ctr_lon, lma_ctr_lat = ds.network_center_longitude.data, ds.network_center_latitude.data

OSW_lat, OSW_lon = 43.454722, -76.506667
Lowville_lat, Lowville_lon = 43.786667, -75.492222
KTYX_lat, KTYX_lon, KTYX_alt = 43.7558, -75.68, 562.0
Smokestack_lat, Smokestack_lon = 43.4592, -76.5312

MET_lat, MET_lon = 43.86273794268999, -75.72739760477387
_144_lat, _144_lon = 43.87550831325867, -75.70463107539224
_083_lat, _083_lon = 43.89425516855627, -75.64199511682814

station_lons = [OSW_lon, Lowville_lon, KTYX_lon, Smokestack_lon]
station_lats = [OSW_lat, Lowville_lat, KTYX_lat, Smokestack_lat]
station_labels = ['Oswego', 'Lowville', 'KTYX', 'Smokestack']
tower_lons = [MET_lon, _144_lon, _083_lon]
tower_lats = [MET_lat, _144_lat, _083_lat]

# Manually set the encoding because someone didn't use UTF-8 like they should have.
#turbines = pd.read_csv('/Users/BenLa/.spyder-py3/turbine_locations.csv', encoding = "ISO-8859-1")
#nearby_turbines = ((turbines.xlong > lma_ctr_lon - 3.0) & (turbines.xlong < lma_ctr_lon + 3.0) & (turbines.ylat > lma_ctr_lat - 3.0) & (turbines.ylat < lma_ctr_lat + 3.0))
#turbines = turbines[nearby_turbines]

DOF = '/Users/BenLa/.spyder-py3/NewYork_Obsticles.Dat'
#DOF = '/Users/BenLa/.spyder-py3/SouthDakota_Obsticles.Dat'

#Bee's\n",
#identify different columns within NewYork_Obsticles.Dat and set them to the corresponding variables
#also skip the first 9000-ish rows as those are all downstate and not apart of our area of interest
specs = [(35, 37), (37, 41), (41, 46), (48, 52), (52, 55), (55, 60), (82, 88), (62, 74)]
StateObs2 = pd.read_fwf(DOF, colspecs = specs, skiprows=4, names=('ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec', 'Elevation (ft)', 'Obsticle'))
#NYObs2 = pd.read_fwf('D:/NewYork_Obsticles.Dat', colspecs = specs, skiprows=9209, skipfooter =1,  names=('ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec', 'Elevation (ft)', 'Obsticle'))

#Bee's
#convert the deg/min/sec for lat and long to decimal format
StateObs2['ylat']=StateObs2['ylatdeg']+StateObs2['ylatmin']/60+StateObs2['ylatsec']/3600
StateObs2['xlong']=(StateObs2['xlongdeg']+StateObs2['xlongmin']/60+StateObs2['xlongsec']/3600)*-1

#Bee's
#get rid of the deg/min/sec columns and send it to a csv file
StateObs2.drop(columns=['ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec']).to_csv(r'/Users/BenLa/.spyder-py3/State_Obsticles.csv')
StateObs3 = pd.read_csv('/Users/BenLa/.spyder-py3/State_Obsticles.csv')


#%%
class AnnotatedLMAPlot(InteractiveLMAPlot):
    # @output.capture()
    def make_plot(self):
        # Use the built-in LMA functionality
        super(AnnotatedLMAPlot, self).make_plot()
        
        # Add our own title
        tlim = self.bounds['t']
        tlim_sub = pd.to_datetime(tlim[0]), pd.to_datetime(tlim[1])
        title = tlim_sub[0].strftime('%Y%m%d %H%M%S') + ' to ' + tlim_sub[1].strftime('%Y%m%d %H%M%S') + ' UTC'
        self.lma_plot.ax_th.set_title(title)
        
        # Add the station positions we defined above, and label them. We only need to add them to the plan-view plot.
        # The other axes are self.lma_plot.ax_th, .ax_lon, and .ax_lat.
        new_artists = []
        art = self.lma_plot.ax_plan.scatter(station_lons, station_lats, s = 1, color='k', zorder=-10, mouseover = True)
        new_artists.append(art)
        for slon,slat,slabel in zip(station_lons, station_lats, station_labels):
            txt_art = self.lma_plot.ax_plan.text(slon, slat, slabel)
            new_artists.append(txt_art)
        
        # Add the wind turbines. We don't worry about filtering to just those in the plot,
        # though we could do that here using 
        #turbine_color=(0,0,1.0,0.5)
        #art = self.lma_plot.ax_plan.scatter(turbines.xlong, turbines.ylat, color=turbine_color, marker='1', zorder=-10)
        #new_artists.append(art)
        #self.data_artists.extend(new_artists)

        #Bee's
        #Add all the other obstacles and turbines from the NewYork_Obsticles File
        obs_colour = (0.0,1.0,0.2,0.3)
        tower_colour = (0.8,0.1,0.3,0.8)
        art = self.lma_plot.ax_plan.scatter(StateObs3.xlong, StateObs3.ylat, color = obs_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(tower_lons, tower_lats, color = tower_colour, marker = '1', zorder=-10)
        new_artists.append(art)
        
class SuperLMAPlot(InteractiveLMAPlot):
    def make_super_plot(self):
        super(SuperLMAPlot, self).make_super_plot()
        # Add our own title
        tlim = self.bounds['t']
        tlim_sub = pd.to_datetime(tlim[0]), pd.to_datetime(tlim[1])
        title = tlim_sub[0].strftime('%Y%m%d %H%M%S') + ' to ' + tlim_sub[1].strftime('%Y%m%d %H%M%S') + ' UTC'
        self.lma_plot.ax_th.set_title(title)
        
        # Add the station positions we defined above, and label them. We only need to add them to the plan-view plot.
        # The other axes are self.lma_plot.ax_th, .ax_lon, and .ax_lat.
        new_artists = []
        art = self.lma_plot.ax_plan.scatter(station_lons, station_lats, s = 1, color='k', zorder=-10, mouseover = True)
        new_artists.append(art)
        for slon,slat,slabel in zip(station_lons, station_lats, station_labels):
            txt_art = self.lma_plot.ax_plan.text(slon, slat, slabel)
            new_artists.append(txt_art)
            
        #Bee's
        #Add all the other obstacles and turbines from the NewYork_Obsticles File
        windmill_array = []
        MET_array = []
        tower_array = []
        bldg_array = []
        utility_array = []
        solar_array = []
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'WINDMILL']:
                #print(StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line])
                windmill_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
                #print(dummy_array[0][0])
        windmill_obs = pd.DataFrame(windmill_array, columns=['xlong', 'ylat', 'obstacles'])
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'MET']:
                MET_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
        MET_obs = pd.DataFrame(MET_array, columns=['xlong', 'ylat', 'obstacles'])
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'TOWER']:
                tower_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
        tower_obs = pd.DataFrame(tower_array, columns=['xlong', 'ylat', 'obstacles'])
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'BLDG']:
                bldg_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
        bldg_obs = pd.DataFrame(bldg_array, columns=['xlong', 'ylat', 'obstacles'])
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'UTILITY POLE']:
                utility_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
        utility_obs = pd.DataFrame(utility_array, columns=['xlong', 'ylat', 'obstacles'])
        for line in StateObs3['Obsticle'].index[StateObs3['Obsticle'] == 'SOLAR PANELS']:
                solar_array.append([StateObs3['xlong'][line], StateObs3['ylat'][line], StateObs3['Obsticle'][line]])
        solar_obs = pd.DataFrame(solar_array, columns=['xlong', 'ylat', 'obstacles'])
        #print(obs)
        turbine_colour = (0.0,1.0,0.2,0.3) #green
        MET_colour = (0.54, 0.5, 0.5, 0.9) #grey
        tower_colour = (0.5,0.1,0.3,0.8) #dark red
        tower_colour2 = (0.1, 0.9, 0.9, 0.8) #light blue
        bldg_colour = (0.3, 0.1, 0.9, 0.9) #purple
        utility_colour = (0.9, 0.7, 0.1, 0.8) #orange
        solar_colour = (0, 0, 0, 1) #black
        art = self.lma_plot.ax_plan.scatter(windmill_obs.xlong, windmill_obs.ylat, color = turbine_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(MET_obs.xlong, MET_obs.ylat, color = MET_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(tower_obs.xlong, tower_obs.ylat, color = tower_colour2, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(bldg_obs.xlong, bldg_obs.ylat, color = bldg_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(tower_lons, tower_lats, color = tower_colour, marker = 'd', zorder=20)
        art = self.lma_plot.ax_plan.scatter(bldg_obs.xlong, bldg_obs.ylat, color = bldg_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(utility_obs.xlong, utility_obs.ylat, color = utility_colour, marker = '1', zorder=-10)
        art = self.lma_plot.ax_plan.scatter(solar_obs.xlong, solar_obs.ylat, color = solar_colour, marker = '1', zorder=-10)
        new_artists.append(art)

radar_data, start_times = lmarad.ReadRadar(["C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_080103_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_080637_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_081311_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_081833_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_082354_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_082916_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_083451_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_084043_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_084649_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_085240_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_085830_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_090420_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_091011_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_091617_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_092223_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_092826_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_093432_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_094037_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_094643_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_095249_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_095855_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_100500_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_101106_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_101710_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_102316_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_102920_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_103526_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_104132_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_104738_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_105344_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_105949_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_110553_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_111159_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_111805_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_112411_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_113017_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_113621_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_114226_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_114831_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_115435_V06"])

gc.collect()
#xsecLatLonAlt = [-104.4, 43.82, -102.9, 44.16, 15] #start lon, start lat, end lon, end lat, maxz
xsecLatLonAlt = [-76.26, 43.64, -75.43, 43.97, 15]

#%% LMA Plot
tlim = pd.to_datetime(starttime).to_pydatetime()-pd.Timedelta(hours=1), pd.to_datetime(endtime).to_pydatetime()-pd.Timedelta(hours=1)

#interactive_lma = AnnotatedLMAPlot(ds, clon=lma_ctr_lon, clat=lma_ctr_lat, tlim=tlim) #, network_data = entln_data)

interactive_lma = SuperLMAPlot(ds, clon=lma_ctr_lon, clat=lma_ctr_lat,tlim=tlim, network_data = entln_data, radar_data = radar_data, points = xsecLatLonAlt)
#%% Flash Statistics <
#Print the horizontal extent of the flash

point_1 = np.where(interactive_lma.this_lma_lon == interactive_lma.this_lma_lon.min())[0]
point_2 = np.where(interactive_lma.this_lma_lon == interactive_lma.this_lma_lon.max())[0]

point_3 = np.where(interactive_lma.this_lma_lat == interactive_lma.this_lma_lat.min())[0]
point_4 = np.where(interactive_lma.this_lma_lat == interactive_lma.this_lma_lat.max())[0]

xy = [[a,b] for (a,b) in zip(interactive_lma.this_lma_lon, interactive_lma.this_lma_lat)]

print("Horizontal Flash Area: ")
points = ConvexHull(xy)
area = 1
for k in range(len(points.vertices)-1): 
    #print(k)
    verticies = points.vertices
    np.append(verticies, points.vertices[0])
    area *= haversine(interactive_lma.this_lma_lat[points.vertices[k]], interactive_lma.this_lma_lon[points.vertices[k]], interactive_lma.this_lma_lat[points.vertices[k+1]], interactive_lma.this_lma_lon[points.vertices[k+1]])/1e3
    #print(area)
print(area**0.5)
#print(points.volume)

#print(haversine(interactive_lma.this_lma_lat[point_1], interactive_lma.this_lma_lon[point_1], interactive_lma.this_lma_lat[point_2], interactive_lma.this_lma_lon[point_2]) * (haversine(interactive_lma.this_lma_lat[point_3], interactive_lma.this_lma_lon[point_3], interactive_lma.this_lma_lat[point_4], interactive_lma.this_lma_lon[point_4]))/1000)

print("Vertical Flash Extent: ")
print(interactive_lma.this_lma_alt.max()-interactive_lma.this_lma_alt.min())

print("Flash Duration: ")
print((float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[0]].strftime('%S.%f'))-float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[-1]].strftime('%S.%f')))*-1)
#print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[0]])
#print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[-1]])
print("Mean Altitude: ")
print(np.mean((interactive_lma.this_lma_alt)))
# Simpler version with no overlays
#first_time = datetime(2022,11,17,0,0,0)
#last_time = datetime(2022,11,21,0,0,0)
#tlim = first_time, last_time
#interactive_lma = InteractiveLMAPlot(ds, tlim = tlim)

#%% STOP
#%% LEE Test Flash Setup
test_flash = {'x': (-75.75868687014818, -75.38774274709453), 
                                 #Longitude from a to b
                             'y': (43.7226927236846, 43.9030922549077),
                                 #Latitude from c to d
                             'z': (0.2372100067138672, 9.768487885779336), 
                                 #Altidue from e to f
                             't': (dt.datetime(2022, 11, 20, 10, 40, 56, 130580), 
                                   dt.datetime(2022, 11, 20, 10, 40, 56, 906314))
                                 #Time from year1, month1, day1, hour1, min1, sec1, microsecond1 to 
                                 #year2, month2, day2, hour2, min2, sec2, microsecond2
                            }

first=0
last=-1

interactive_lma.lma_plot.ax_th.set_xlim(test_flash['t'], emit=False)
interactive_lma.lma_plot.ax_th.set_ylim(test_flash['z'], emit=False)
interactive_lma.lma_plot.ax_plan.set_xlim(test_flash['x'], emit=False)
interactive_lma.lma_plot.ax_plan.set_ylim(test_flash['y'])

#%% Zoom the box to fit FED plot
test_flash = {'x': (-77.1605462881235, -75.41592407226562),
              'y': (43.38042449951172, 44.09492978182706),
              'z': (0.2372100067138672, 10.013169921875),
              't': (dt.datetime(2022, 11, 17, 7, 0),
                    dt.datetime(2022, 11, 21, 11, 0))}

first=0
last=-1

interactive_lma.lma_plot.ax_th.set_xlim(test_flash['t'], emit=False)
interactive_lma.lma_plot.ax_th.set_ylim(test_flash['z'], emit=False)
interactive_lma.lma_plot.ax_plan.set_xlim(test_flash['x'], emit=False)
interactive_lma.lma_plot.ax_plan.set_ylim(test_flash['y'])

#%% Velocity Plot <
distance_from_origin = []
time_from_origin = []
altitude = []
i=0
m=-10

first=0
last=-1

while (i < (interactive_lma.this_lma_lat.size)):
    distance_from_origin.append(haversine(interactive_lma.this_lma_lat[first],interactive_lma.this_lma_lon[first],interactive_lma.this_lma_lat[i],interactive_lma.this_lma_lon[i]))
    time_from_origin.append((float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[first]].strftime('%S.%f'))-float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[i]].strftime('%S.%f')))*-1)
    altitude.append(interactive_lma.this_lma_alt[i])
    i+=1

fig = plt.figure(figsize=(10, 10))
while (m < 10):
    x = np.linspace(0, (float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[first]].strftime('%S.%f'))-float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[last]].strftime('%S.%f')))*-1+10, 100)
    y = x * 2 * 10**4
    plt.plot(x+(m/10), y, color = 'b')
    yy = x * 10**5
    plt.plot(x+(m/10), yy, color = 'r')
    yyy = x * 10**6
    plt.plot(x+(m/10), yyy, color = 'g')
    m+=1

x = np.linspace(0, (float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[first]].strftime('%S.%f'))-float(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[last]].strftime('%S.%f')))*-1+10, 100)
y = x * 2 * 10**4
plt.plot(x, y, color = 'b', label = 'Positive Leader')
yy = x * 10**5
plt.plot(x, yy, color = 'r', label = 'Negative Leader')
yyy = x * 10**6
plt.plot(x, yyy, color = 'g', label = 'Dart Leader')
    
plt.ylim(0, max(distance_from_origin)+1000)
plt.xlim(-0.05, max(time_from_origin)+0.1) 

sc = plt.scatter(time_from_origin, distance_from_origin, zorder=10, c=altitude, cmap="cool", vmin=0, vmax=15)
plt.colorbar(sc, label='Altitude (km)', spacing='proportional', ticks=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
plt.legend(title = "Leader Type Key", loc = 2, fontsize=12, title_fontsize=16, framealpha=1)
plt.title("Lightning Leader Speed", size=24)
plt.xlabel("Time from Origin (s)", size=12)
plt.ylabel("Distance from Origin (m)", size=12)

plt.show()

#%% Radar Download
# Set the WSR-88D radar name
radar_id = 'KTYX'
#radar_id = 'KTLX'
#radar_id = 'KUDX'
begin = dt.datetime(2022,11,20,11)
spacing = [0, 10]

## Set the plotting criteria!
min_stations = 6 # more stations = more confident it's a good solution
max_chi = 1 # lower reduced chi^2 = more confident it's a good solution
tbuffer = 6*60 # How many seconds AFTER the radar scan time do you want to plot LMA data?
max_dist = 2.5e3 # meters, max distance of an LMA source to a radar cross-section plane to be plotted

min_events_per_flash = 5 # Minimum number of sources per flash

#Create datetime objects for the start & end times to download:
start=begin - dt.timedelta(minutes=spacing[0])
end = start + dt.timedelta(minutes=spacing[1])

# Create a folder to store the data
#downloadloc = f'{start:%Y%m%d%H}_{radar_id}'

# Determine What scans are avaiable for the radar site and times listed
scans = conn.get_avail_scans_in_range(start, end, radar_id)

print("There are {} scans available between {} and {}".format(len(scans), start, end))
print(scans[0:4])

#download the data
results = conn.download(scans, downloadloc)

files = []
files = sorted(glob.glob(f'{downloadloc}/*V06'))

#####
# # OR if already downloaded just replace the above with:
#####
# files = sorted(glob.glob('path/to/88Ddata/*V06'))
#%% Radar Plan View <

# Select a file of interest
which_file = files[0]
file_time = [dt.datetime.strptime(which_file.split('/')[-1],"KUDX%Y%m%d_%H%M%S_V06")]

# Read it into a pyart object and find start and end times to search the LMA data
radar_data =  pyart.io.read(which_file)
start = begin - dt.timedelta(seconds = tbuffer)
end   = file_time[0] + dt.timedelta(seconds = tbuffer)

mx,xx = -104.4,-102.9 # left and right end points (x, km) 
my,xy = 43.8,44.16 # up and down end point (y, km)

#mx, xx = DOW_lon, END_lon[0]
#my, xy = DOW_lat, END_lat[0]

max_z = 15
max_dist_deg = 0.031 # degrees, max distance of an LMA source to a radar grid centroid to be plotted, roughtly 2.5km in midlatitudes
radarx = radar_data.longitude['data'][0]
radary = radar_data.latitude['data'][0]

mx_rad, xx_rad = haversine(my, mx, my, radarx)/1e3, haversine(xy, xx, xy, radarx)/1e3
my_rad, xy_rad = haversine(radary, mx, my, mx)/1e3, haversine(radary, xx, xy, xx)/1e3

if negativeaxis(mx, radar_data): mx_rad = mx_rad*-1
if negativeaxis(my, radar_data): my_rad = my_rad*-1
if negativeaxis(xx, radar_data): xx_rad = xx_rad*-1
if negativeaxis(xy, radar_data): xy_rad = xy_rad*-1

nn = int(((xx_rad-mx_rad)**2+(xy_rad-my_rad)**2)**0.5/0.1) 
des_x = np.linspace(mx_rad,xx_rad,nn)*1e3
des_y = np.linspace(my_rad,xy_rad,nn)*1e3
des_z = np.arange(0,max_z+0.1,0.1)*1e3

desx_grid,desz_grid = np.meshgrid(des_x,des_z)
desy_grid,desz_grid = np.meshgrid(des_y,des_z)
new_x = np.arange(0,np.shape(des_x)[0]*0.1 ,0.1)*1e3

keep_sweeps_pol  = [] # Holding spot for good polarimetric variable (non-NaN ZDR) sweeps
keep_sweeps_trad = [] # Holding spot for good traditional variable (non-NaN velocity) sweeps
for sweep in radar_data.sweep_number['data']:
    if not (np.size(radar_data.extract_sweeps([sweep]).fields['differential_reflectivity']['data'].mask)-
             np.sum(radar_data.extract_sweeps([sweep]).fields['differential_reflectivity']['data'].mask))==0:
        #print ('Keeping polarimetric sweep number: ', sweep)
        keep_sweeps_pol+=[sweep]
    if not (np.size(radar_data.extract_sweeps([sweep]).fields['velocity']['data'].mask)-
             np.sum(radar_data.extract_sweeps([sweep]).fields['velocity']['data'].mask))==0:
        #print ('Keeping traditional sweep number: ', sweep)
        keep_sweeps_trad+=[sweep]

# Keep only the good sweeps for each set of radar variables
radar_pol = radar_data.extract_sweeps(keep_sweeps_pol)
radar_trad = radar_data.extract_sweeps(keep_sweeps_trad)

des_tree_pol = spatial.KDTree(np.array([radar_pol.gate_x['data'].ravel(),radar_pol.gate_y['data'].ravel(), radar_pol.gate_z['data'].ravel()]).T)

des_tree_trad = spatial.KDTree(np.array([radar_trad.gate_x['data'].ravel(),radar_trad.gate_y['data'].ravel(), radar_trad.gate_z['data'].ravel()]).T)

# Find the radar gate closest to the cross-section grid
dists_pol, indext_pol = des_tree_pol.query(np.array([desx_grid.ravel(), desy_grid.ravel(), desz_grid.ravel()]).T)

# Find the radar gate closest to the cross-section grid
dists_trad, indext_trad = des_tree_trad.query(np.array([desx_grid.ravel(),  desy_grid.ravel(),  desz_grid.ravel()]).T)

# Get the x,y,z locations of the radar gates for later
# Theoretically these should match, but let's not make any assumptions
rx_pol, ry_pol, rz_pol  = radar_pol.get_gate_x_y_z(0)
rx_trad,ry_trad,rz_trad = radar_trad.get_gate_x_y_z(0)

sweep=1 # sweep 0 = lowest level
fig = plt.figure(figsize=(15,10))
ax = plt.subplot(projection = ccrs.PlateCarree())
display = pyart.graph.RadarMapDisplay(radar_data)
display.plot_ppi_map('reflectivity', sweep, vmin = -20, vmax=60, cmap='pyart_HomeyerRainbow', min_lat=radar_data.latitude['data'][0]-2, max_lat=radar_data.latitude['data'][0]+2, min_lon=radar_data.longitude['data'][0]-2, max_lon=radar_data.longitude['data'][0]+2, ax = ax)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
ax.add_feature(cfeature.BORDERS)
                     
plt.plot([mx, xx], [my, xy], color = '#6a8f9c', linewidth = 2)

min_time = min(set(interactive_lma.this_lma_time.values))
max_time = max(set(interactive_lma.this_lma_time.values))

fruit = 0
veggies = 0
plop = []
save_xsec_lma = []
trash_xsec_lma = []

while fruit < len(interactive_lma.this_lma_alt):
        if interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') >= start.strftime('%H:%M:%S.%f') and interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') <= end.strftime('%H:%M:%S.%f'):
                #print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f'))
                plop.append(np.cross((xx, xy)-np.array([mx, my]), np.array([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit]])-np.array([mx, my]))/norm((xx, xy)-np.array([mx, my])))
                #print(cake)
                #print(d)
                if abs(plop[veggies]) < max_dist_deg:
                        save_xsec_lma.append([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], interactive_lma.this_lma_alt[fruit]])
                        plt.scatter(interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], color = 'b', s = 5)
                else:
                        trash_xsec_lma.append([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], interactive_lma.this_lma_alt[fruit]])
                        plt.scatter(interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], color = 'k', s = 5)
                veggies += 1
        fruit += 1

included_values = [m for m in entln_data['datetime'].values if min_time < m < max_time]
indexed_values = np.where(entln_data['datetime'].isin(included_values))[0]

Positive_IC_lat, Negative_IC_lat, Positive_CG_lat, Negative_CG_lat, Positive_IC_lon, Negative_IC_lon, Positive_CG_lon, Negative_CG_lon = [], [], [], [], [], [], [], []

"""
Positive_IC_lat = []
Negative_IC_lat = []
Positive_CG_lat = []
Negative_CG_lat = []
Positive_IC_lon = []
Negative_IC_lon = []
Positive_CG_lon = []
Negative_CG_lon = []
"""
Positive_IC_alt = [max_z-1] * indexed_values.size
Negative_IC_alt = [max_z-1] * indexed_values.size
Positive_CG_alt = [5] * indexed_values.size
Negative_CG_alt = [5] * indexed_values.size

# # Make sure it looks reasonable from top-down view
#plt.scatter((interactive_lma.this_lma_lon), (interactive_lma.this_lma_lat), color='k',s=5) # all sources meeting criteria
j = 0
while j < indexed_values.size:
    if (entln_data['ic height'][indexed_values].values[j] == 8000):
        #print('IC Flash')
        #IC Flash
        if entln_data['peak_current_kA'][indexed_values].values[j] > 0:
            #print('Positive IC')
            Positive_IC_lat.append(entln_data['latitude'][indexed_values].values[j])
            Positive_IC_lon.append(entln_data['longitude'][indexed_values].values[j])
            plt.scatter(entln_data['longitude'][indexed_values].values[j], entln_data['latitude'][indexed_values].values[j], marker = '^', color='blue',s=5)
        else:
            #print('Negative IC')
            Negative_IC_lat.append(entln_data['latitude'][indexed_values].values[j])
            Negative_IC_lon.append(entln_data['longitude'][indexed_values].values[j])
            plt.scatter(entln_data['longitude'][indexed_values].values[j], entln_data['latitude'][indexed_values].values[j], marker = '^', color='red',s=5)
    else:
        #print('CG Flash')
        #CG Flash
        if entln_data['peak_current_kA'][indexed_values].values[j] > 0:
            #print('Positive CG')
            Positive_CG_lat.append(entln_data['latitude'][indexed_values].values[j])
            Positive_CG_lon.append(entln_data['longitude'][indexed_values].values[j])
            plt.scatter(entln_data['longitude'][indexed_values].values[j], entln_data['latitude'][indexed_values].values[j], marker = 'v', color='blue',s=5)
        else:
            #print('Negative CG')
            Negative_CG_lat.append(entln_data['latitude'][indexed_values].values[j])
            Negative_CG_lon.append(entln_data['longitude'][indexed_values].values[j])
            plt.scatter(entln_data['longitude'][indexed_values].values[j], entln_data['latitude'][indexed_values].values[j], marker = 'v', color='red',s=5)
    j+=1
plt.show()

#%% Specific Differential Phase (Only Run Once)
kdpdata = pyart.retrieve.kdp_vulpiani(radar_trad,  phidp_field ='differential_phase',  band ='S' , windsize = 34)
kdpdata = kdpdata[0]['data'] #Retreive the kdp 'data' from the the kdpdata that was calculated using the Vulpiani method 
mask = np.logical_and(kdpdata > -0.01, kdpdata < 0.01) #mask the values from -0.01 to 0.01 so they do not plot
kdpdata = np.where(mask, np.nan, kdpdata) #Apply the mask to kdpdata
radar_trad.add_field('specific_differential_phase_hv', {'data': kdpdata.data})  #Add a new field to the radar data dictionary with the kdp data
#%% Radar Cross Sections <
# Look for LMA sources near the radar grid points

lma_alt = interactive_lma.this_lma_alt
new_lma_x = interactive_lma.this_lma_lon
new_lma_y = interactive_lma.this_lma_lat

t = 0
o = 0
m = 0
a = 0
s = 0

new_lma_r = []
new_lma_alt = []
new_netwPIC_r = []
new_netwNIC_r = []
new_netwPCG_r = []
new_netwNCG_r = []

while t < len(save_xsec_lma)-1:
    new_lma_r.append(haversine(my, mx, save_xsec_lma[t][1], save_xsec_lma[t][0])/1e3)
    new_lma_alt.append(save_xsec_lma[t][2])
    #new_lma_r.append(haversine(my, mx, new_lma_y[t], new_lma_x[t])/1e3)
    t+=1
while o <= len(Positive_IC_lat)-1:
    if len(Positive_IC_lat) > 0: new_netwPIC_r.append(haversine(my, mx, Positive_IC_lat[o], Positive_IC_lon[o])/1e3)
    o+=1
while m <= len(Negative_IC_lat)-1:
    if len(Negative_IC_lat) > 0: new_netwNIC_r.append(haversine(my, mx, Negative_IC_lat[m], Negative_IC_lon[m])/1e3)
    m+=1
while a <= len(Positive_CG_lat)-1:
    if len(Positive_CG_lat) > 0: new_netwPCG_r.append(haversine(my, mx, Positive_CG_lat[a], Positive_CG_lon[a])/1e3)
    a+=1
while s <= len(Negative_CG_lat)-1:
    if len(Negative_CG_lat) > 0: new_netwNCG_r.append(haversine(my, mx, Negative_CG_lat[s], Negative_CG_lon[s])/1e3)
    s+=1

fig = plt.figure(figsize=(28,16))
#Reflectivity
ax1 = fig.add_subplot(321)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['reflectivity']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(-20,72,1),cmap='pyart_HomeyerRainbow')
plt.colorbar(label='Reflectivity (dBZ)')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
plt.ylim(0, max_z)
#Velocity
ax2 = fig.add_subplot(322)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['velocity']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(-40,40,1),cmap='NWSVel')
plt.colorbar(label='Velocity (m/s)')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
plt.ylim(0, max_z)
#Spectrum Width
ax3 = fig.add_subplot(323)
colors1 = plt.cm.binary_r(np.linspace(0.2,0.8,33))
colors2 = plt.cm.gnuplot_r(np.linspace(0.,0.7,100))
colors = np.vstack((colors1, colors2[10:121]))
swcolours = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['spectrum_width']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(0,14,0.1),cmap='NWS_SPW')
plt.colorbar(label='Spectrum Width')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
#Differential Reflectivity
ax4 = fig.add_subplot(324)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_pol.fields['differential_reflectivity']['data'].ravel()[indext_pol].reshape(np.shape(desz_grid)),levels = np.arange(-4,8,0.1),cmap='ChaseSpectral')
plt.colorbar(label='Differential Reflectivity (ZDR)')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
plt.ylim(0, max_z)
#Correlation Coefficient
ax5 = fig.add_subplot(325)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_pol.fields['cross_correlation_ratio']['data'].ravel()[indext_pol].reshape(np.shape(desz_grid)),levels = np.arange(0,1.1,0.01),cmap='SCook18')
plt.colorbar(label='Correlation Coefficient')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
plt.ylim(0, max_z)
#Specific Differential Phase
ax6 = fig.add_subplot(326)
plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, kdpdata.ravel()[indext_pol].reshape(np.shape(desz_grid)),np.arange(-1,4.1,0.1),cmap=swcolours)
plt.colorbar(label='Specific Differential Phase (KDP)')
plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
if len(new_netwPIC_r) > 0: plt.scatter(new_netwPIC_r, [max_z-10]*len(new_netwPIC_r), color='b', marker='^', s=25)
if len(new_netwNIC_r) > 0: plt.scatter(new_netwNIC_r, [max_z-10]*len(new_netwNIC_r), color='r', marker='^', s=25)
if len(new_netwPCG_r) > 0: plt.scatter(new_netwPCG_r, [1]*len(new_netwPCG_r), color='b', marker='v', s=25)
if len(new_netwNCG_r) > 0: plt.scatter(new_netwNCG_r, [1]*len(new_netwNCG_r), color='r', marker='v', s=25)
plt.xlim(0, np.max(new_x)/1e3)
plt.ylim(0, max_z)
         
         
plt.ylabel('Altitude (km AGL)')
plt.xlabel('Distance (km)')
plt.tight_layout()

#%% Radar Plan View with DOW
min_stations = 6 # more stations = more confident it's a good solution
max_chi = 1 # lower reduced chi^2 = more confident it's a good solution
tbuffer = 2.5*60 # How many seconds before and after the radar scan time do you want to plot LMA data?
max_dist_deg = 0.031 # degrees, max distance of an LMA source to a radar grid centroid to be plotted, roughtly 2.5km in midlatitudes

# Set max range for DOW RHI
rng_rhi = 55

# Start with the DOW file of interest
dow_file = 'C:/Users/BenLa/Downloads/cfrad.20221120_110610.994_DOW7_v135_s00_az359.90_RHI.nc'
# Read into pyart
dow_pyart = pyart.io.read(dow_file)
gatefilter = pyart.filters.GateFilter(dow_pyart)
gatefilter.exclude_below('DBZHCC', 10)
# Find x,y,z of gates
rx,ry,rz = dow_pyart.get_gate_x_y_z(0)
# Find the starttime
dow_starttime = dt.datetime.strptime(dow_pyart.time['units'].split(' ')[-1],
                     '%Y-%m-%dT%H:%M:%SZ')+dt.timedelta(seconds = dow_pyart.time['data'][0])

ref_field = ma.getdata(dow_pyart.fields['DBZHCC']['data'])
ref_gt_0 = np.ma.where(ref_field > 10, ref_field, -40)
vel_gt_0 = ma.masked_array(dow_pyart.fields['VEL']['data'], mask = ref_field < 10)
sw_gt_0 = ma.masked_array(dow_pyart.fields['WIDTH']['data'], mask = ref_field < 10)
zdr_gt_0 = ma.masked_array(dow_pyart.fields['ZDRC']['data'], mask = ref_field < 10)
cc_gt_0 = ma.masked_array(dow_pyart.fields['RHOHV']['data'], mask = ref_field < 10)
kdp_gt_0 = ma.masked_array(dow_pyart.fields['KDP']['data'], mask = ref_field < 10)
#print(ref_gt_0)
mask_dict_ref = {"data": ref_gt_0,"units": "dBZ","long_name": "reflectivity_mask","_FillValue": ref_gt_0.fill_value,"standard_name": "reflectivity_mask"}
mask_dict_vel = {'data': vel_gt_0, 'Long_name': 'velocity_mask', '_FillValue': ref_gt_0.fill_value, 'standard_name': 'velocity_mask'}
mask_dict_sw = {'data': sw_gt_0, 'Long_name': 'spectrum_width_mask', 'standard_name': 'spectrum_width_mask'}
mask_dict_zdr = {'data': zdr_gt_0, 'Long_name': 'differential_reflectivity_mask', 'standard_name': 'differential_reflectivity_mask'}
mask_dict_cc = {'data': cc_gt_0, 'Long_name': 'correlation_coefficient_mask', 'standard_name': 'correlation_coefficient_mask'}
mask_dict_kdp = {'data': kdp_gt_0, 'Long_name': 'specific_differential_phase_mask', 'standard_name': 'specific_differential_phase_mask'}
dow_pyart.add_field("reflectivity_mask", mask_dict_ref, replace_existing=True)
dow_pyart.add_field('velocity_mask', mask_dict_vel, replace_existing=True)
dow_pyart.add_field('spectrum_width_mask', mask_dict_sw, replace_existing=True)
dow_pyart.add_field('differential_reflectivity_mask', mask_dict_zdr, replace_existing=True)
dow_pyart.add_field('correlation_coefficient_mask', mask_dict_cc, replace_existing=True)
dow_pyart.add_field('specific_differential_phase_mask', mask_dict_kdp, replace_existing=True)


# Pull the azimuth name from the file name (degrees)
az_name = dow_file.split('/')[-1].split('az')[1].split('_')[0]
# And also the mode reported in the file - convert this to form used later
# for calculations (in radians)
az_rad  = np.deg2rad(stats.mode( dow_pyart.azimuth["data"])[0]-90)

print (dow_starttime, " at azimuth= ", az_name)

# Find the start and end times for pulling LMA data
start = dow_starttime - dt.timedelta(seconds=tbuffer)
end   = dow_starttime + dt.timedelta(seconds=tbuffer)

lma_alt = interactive_lma.this_lma_alt
new_lma_x = interactive_lma.this_lma_lon
new_lma_y = interactive_lma.this_lma_lat

DOW_lon = dow_pyart.get_gate_lat_lon_alt(sweep=0)[1][0,0]
DOW_lat = dow_pyart.get_gate_lat_lon_alt(sweep=0)[0][0,0]
g1_x, g1_y = pyart.core.geographic_to_cartesian_aeqd(DOW_lon,DOW_lat, 
                                                    radar_data.longitude['data'][0], 
                                                    radar_data.latitude['data'][0], 
                                                    R=6370997.0)
g2_x =  np.cos(az_rad)*rng_rhi
g2_y = -np.sin(az_rad)*rng_rhi
# And the lat/lon of that point
END_lon, END_lat = pyart.core.cartesian_to_geographic_aeqd(g1_x+g2_x*1e3, g1_y+g2_y*1e3,
                                                        radar_data.longitude['data'][0], 
                                                        radar_data.latitude['data'][0], 
                                                        R=6370997.0)

mx_DOW, xx_DOW = haversine(DOW_lat, DOW_lon, DOW_lat, radarx)/1e3, haversine(END_lat, END_lon, END_lat, radarx)/1e3
my_DOW, xy_DOW = haversine(radary, DOW_lon, DOW_lat, DOW_lon)/1e3, haversine(radary, END_lon, END_lat, END_lon)/1e3

nn = int(((xx_DOW-mx_DOW)**2+(xy_DOW-my_DOW)**2)**0.5/0.1) 
DOW_des_x = np.linspace(mx_DOW,xx_DOW,nn)*1e3
DOW_des_y = np.linspace(my_DOW,xy_DOW,nn)*1e3
DOW_des_z = np.arange(0,max_z+0.1,0.1)*1e3

desx_grid,desz_grid = np.meshgrid(DOW_des_x,DOW_des_z)
desy_grid,desz_grid = np.meshgrid(DOW_des_y,DOW_des_z)
DOWnew_x = np.arange(0,np.shape(des_x)[0]*0.1 ,0.1)*1e3

sweep=1 # sweep 0 = lowest level
fig = plt.figure(figsize=(15,10))
ax = plt.subplot(projection = ccrs.PlateCarree())
display = pyart.graph.RadarMapDisplay(radar_data)
display.plot_ppi_map('reflectivity', sweep, vmin = -20, vmax=60, cmap='pyart_HomeyerRainbow', min_lat=radar_data.latitude['data'][0]-2, max_lat=radar_data.latitude['data'][0]+2, min_lon=radar_data.longitude['data'][0]-2, max_lon=radar_data.longitude['data'][0]+2, ax = ax)
ax.set_xlim(-76.28, -75.51)
ax.set_ylim(43.33, 44.07)

ax.arrow(DOW_lon, DOW_lat, 
         END_lon[0]-DOW_lon,END_lat[0]-DOW_lat,
         color='black', linewidth=1, head_width=0.01)
ax.scatter(DOW_lon, DOW_lat, color='k', marker='x')

bread = 0
cake = 0
d = []
save_lma = []
trash_lma = []

while bread < len(interactive_lma.this_lma_alt):
        if interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[bread]].strftime('%H:%M:%S.%f') >= start.strftime('%H:%M:%S.%f') and interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[bread]].strftime('%H:%M:%S.%f') <= end.strftime('%H:%M:%S.%f'):
                #print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[bread]].strftime('%H:%M:%S.%f'))
                d.append(np.cross((END_lon[0], END_lat[0])-np.array([DOW_lon, DOW_lat]),
                           np.array([interactive_lma.this_lma_lon[bread], interactive_lma.this_lma_lat[bread]])
                           -np.array([DOW_lon, DOW_lat]))/norm((END_lon[0], END_lat[0])-np.array([DOW_lon, DOW_lat])))
                #print(cake)
                #print(d)
                if abs(d[cake]) < max_dist_deg:
                        save_lma.append([interactive_lma.this_lma_lon[bread], interactive_lma.this_lma_lat[bread], interactive_lma.this_lma_alt[bread]])
                        plt.scatter(interactive_lma.this_lma_lon[bread], interactive_lma.this_lma_lat[bread], color = 'b', s = 5)
                else:
                        trash_lma.append([interactive_lma.this_lma_lon[bread], interactive_lma.this_lma_lat[bread], interactive_lma.this_lma_alt[bread]])
                        plt.scatter(interactive_lma.this_lma_lon[bread], interactive_lma.this_lma_lat[bread], color = 'k', s = 5)
                cake += 1
        bread += 1

ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
ax.add_feature(cfeature.BORDERS)

fig.text(0.17, 0.89, '(a)', fontsize = 22, fontweight = 'bold')

plt.tight_layout()
plt.show()

#%% DOW Cross Section

wheat = 0
DOW_new_lma_r = []
DOW_lma_alt = []
while wheat < len(save_lma):
    DOW_new_lma_r.append(haversine(DOW_lat, DOW_lon, save_lma[wheat][1], save_lma[wheat][0])/1e3)
    DOW_lma_alt.append(save_lma[wheat][2])
    wheat+=1

fig = plt.figure(figsize = (28, 22))
plt.title('DOW Scan ' + dow_file.split('/')[-1].split('_')[0] + "\n" + dow_file.split('/')[-1].split('_')[1] + ' ' + dow_file.split('/')[-1].split('_')[-2] + "\n", fontsize=30)
plt.axis('off')
ax1 = fig.add_subplot(321)
# The DOW RHI
#display.plot('reflectivity', vmin = -20, vmax = 62, cmap = 'pyart_HomeyerRainbow', gatefilter = gatefilter, ax = ax1)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['reflectivity_mask']['data'],
            cmap = 'pyart_HomeyerRainbow', levels = np.arange(-20,62,2),
            gatefilter = gatefilter)
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label = 'Reflectivity (dBZ)')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label = 'Reflectivity (dBZ)', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.02, 0.88, '(b)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

ax2 = fig.add_subplot(322)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['velocity_mask']['data'],
            cmap = 'pyart_NWSVel', levels = np.arange(-40,40,1))
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label='Velocity (m/s)')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Velocity (m/s)', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.51, 0.88, '(c)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

ax3 = fig.add_subplot(323)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['spectrum_width_mask']['data'],
            cmap = 'NWS_SPW', levels = np.arange(0,14,0.1))
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label='Spectrum Width')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Spectrum Width', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.02, 0.575, '(d)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

ax4 = fig.add_subplot(324)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['differential_reflectivity_mask']['data'],
            cmap = 'ChaseSpectral', levels = np.arange(-4,8,0.1))
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label='Differential Reflectivity (ZDR)')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Differential Reflectivity (ZDR)', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.51, 0.575, '(e)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

ax5 = fig.add_subplot(325)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['correlation_coefficient_mask']['data'],
            cmap = 'SCook18', levels = np.arange(0,1.1,0.01))
ax5.tick_params(axis='x', labelsize=16)
ax5.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label='Correlation Coefficient')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Correlation Coefficient', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.02, 0.27, '(f)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

ax6 = fig.add_subplot(326)
colors1 = plt.cm.binary_r(np.linspace(0.2,0.8,33))
colors2 = plt.cm.gnuplot_r(np.linspace(0.,0.7,100))
colors = np.vstack((colors1, colors2[10:121]))
swcolours = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
plt.contourf((rx**2+ry**2)**0.5/1e3,rz/1e3,
            dow_pyart.fields['specific_differential_phase_mask']['data'],
            cmap = swcolours, levels = np.arange(-1,4.1,0.1))
ax6.tick_params(axis='x', labelsize=16)
ax6.tick_params(axis='y', labelsize=16)
cbar = plt.colorbar(label='Specific Differential Phase (KDP)')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Specific Differential Phase (KDP)', size = 16)
plt.scatter(DOW_new_lma_r, DOW_lma_alt, color='k', s=10)
fig.text(0.51, 0.27, '(g)', fontsize = 22, fontweight = 'bold')
plt.ylim(0, max_z)
plt.xlim(0, rng_rhi)

plt.tight_layout()

#%% Create the Averaged Radar Grid to Plot With

radar_pickle = files[0]+str(len(files))+'.pickle'
reflectivity_grids = []

if os.path.exists(radar_pickle):
    with open(radar_pickle, 'rb') as pickle_file:
        grided_radar = pickle.load(pickle_file) 
else:
        radar_list, start_times = ReadRadar(files)
        grided_radar = GridRadar(radar_list)

        with open(radar_pickle, 'wb') as pickle_file:
                pickle.dump(grided_radar, pickle_file)

filler_grid = np.full_like(grided_radar[0].fields['reflectivity']['data'].data, -10)
for number in grided_radar[:-1]:
        #filler_grid = np.full_like(number.fields["reflectivity"]['data'].data, -10)
        #print(height, row, column)
        height = 0
        while height < grided_radar[0].fields['reflectivity']['data'].shape[0]:
                #print('height')
                #print(height)
                row = 0
                while row < grided_radar[0].fields['reflectivity']['data'].shape[1]:
                        #print('row')
                        #print(row)
                        column = 0
                        while column < grided_radar[0].fields['reflectivity']['data'].shape[2]:
                                #print("column")
                                #print(column)
                                if filler_grid[0][row][column] < number.fields['reflectivity']['data'].data[height][row][column]:
                                        filler_grid[0][row][column] = number.fields['reflectivity']['data'].data[height][row][column]
                                        #print(filler_grid[row][column])
                                column += 1
                        row += 1
                height += 1
        filler_grid[0] = np.select([filler_grid[0] < 1], [np.nan], filler_grid[0])
        composite_dict = {"data": filler_grid, 'units': 'dBZ', 'long_name': 'composite_reflectivity', '_FillValue': filler_grid, 'standard_name': 'composite_reflectivity'}
        number.add_field('composite_reflectivity', composite_dict, replace_existing=True)
        
        reflectivity_grids.append(10**(number.fields["reflectivity"]['data'][0]/10))
        
        #reflectivity_grids.append((number.fields["reflectivity"]['data'][0]))
        #print(number.fields['reflectivity']['data'][0])

save = np.full_like(reflectivity_grids[0], 0)
divisor = np.full_like(reflectivity_grids[0], 0)
for amount in reflectivity_grids:
        for cols in range(amount.shape[1]):
                for rows in range(amount.shape[0]):
                        if amount[cols][rows]:
                                save[cols][rows] += amount[cols][rows]
                         
save = 10 * np.log10(save/len(reflectivity_grids))

fig = plt.figure(figsize=(20, 20))
ax = plt.subplot(projection = ccrs.PlateCarree())                                                    #left, right, bottom, top
ax.imshow(save,vmin = -20, vmax=60, cmap='pyart_HomeyerRainbow',extent=(grided_radar[0].origin_longitude['data'][0]-1.5, grided_radar[0].origin_longitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]-1.5))
ax.invert_yaxis() 
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES, edgecolor = 'blue')

plt.show()

#%% Create a Radar Frequency Grid to Plot With

files = ["C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_173520_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_174026_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_174532_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_175024_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_175516_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_180021_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_180543_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_181049_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_181625_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_182147_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_182721_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_183326_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_183933_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_184538_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_185143_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_185748_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_190354_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_190945_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_191507_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_192027_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_192549_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_193055_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_193546_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_194037_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_194514_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_194950_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_195426_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_195904_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_200340_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_200817_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_201254_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_201730_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_202207_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_202643_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_203120_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_203552_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_204027_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_204503_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_205010_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_205447_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_205924_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_210400_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_210837_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_211314_V06",
"C:/Users/BenLa/.spyder-py3/2022111808_KTYX/KTYX20221120_211750_V06"]

radar_pickle = files[0]+str(len(files))+'.pickle'
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
frequency_grids = []

if os.path.exists(radar_pickle):
    with open(radar_pickle, 'rb') as pickle_file:
        radar_frequency_files = pickle.load(pickle_file) 
else:
        radar_frequency_files, radar_frequency_start = lmarad.ReadRadar(files)

        with open(radar_pickle, 'wb') as pickle_file:
                pickle.dump(radar_frequency_files, pickle_file)


for number in radar_frequency_files[:-1]:
        #print(number)
        frequency_grids.append((number.fields["reflectivity"]['data'][0]))
        
        
fsave = np.full_like(frequency_grids[0], 0)
for famount in frequency_grids:
        for fcols in range(famount.shape[1]):
                for frows in range(famount.shape[0]):
                        if famount[fcols][frows] > 20:
                                fsave[fcols][frows] += 1
fsave = fsave/np.max(fsave)
                                
fig = plt.figure(figsize=(20, 20))
ax = plt.subplot(projection = ccrs.PlateCarree())
im = ax.imshow(fsave, vmin = 0, vmax = np.max(fsave), cmap = 'HomeyerRainbow', alpha = 0.5, extent=(radar_frequency_files[0].origin_longitude['data'][0]-1.5, radar_frequency_files[0].origin_longitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]-1.5))                                
ok = ax.contour(fsave, vmin = 0, vmax = np.max(fsave)/np.max(fsave), levels = levels, cmap = 'HomeyerRainbow', origin = 'upper', extent=(radar_frequency_files[0].origin_longitude['data'][0]-1.5, radar_frequency_files[0].origin_longitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]-1.5))
ax.invert_yaxis()
fig.colorbar(im)
fig.colorbar(ok)
ax.add_feature(COUNTIES, facecolor = 'none', edgecolor = 'gray')
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES, edgecolor = 'blue')

plt.show()


#%% Flash Extent Density Plot

grid_size = 0.025 #0.025 degress lat/lon
xedge = np.arange(int(interactive_lma.bounds['x'][0])-1, int(interactive_lma.bounds['x'][1])+1, grid_size)
yedge = np.arange(int(interactive_lma.bounds['y'][0])-1, int(interactive_lma.bounds['y'][1])+1, grid_size)


eventx = interactive_lma.this_lma_lon
eventy = interactive_lma.this_lma_lat
eventt = interactive_lma.this_lma_time

event_center = (np.mean([interactive_lma.bounds['x'][0], interactive_lma.bounds['x'][1]]),
                np.mean([interactive_lma.bounds['y'][0], interactive_lma.bounds['y'][1]]))

#grid_edge = [event_center[0]-zoom, event_center[0]+zoom, event_center[1]+zoom, event_center[1]-zoom] #left, right, top, bottom
#gridx = np.linspace(start = (event_center[0]-zoom), stop = (event_center[0]+zoom), num = grid_number)
#gridy = np.linspace(start = (event_center[1]+zoom), stop = (event_center[1]-zoom), num = grid_number)
#mesh_gridx, mesh_gridy = np.meshgrid(gridx, gridy)

grid_points, xedges, yedges = np.histogram2d(eventy, eventx, bins = (yedge, xedge))

#grided_data = griddata((eventx, eventy), grid_points, (mesh_gridx, mesh_gridy), method = 'linear')

grid_fig = plt.figure(figsize=(10, 4.5))
ax = plt.subplot(projection = ccrs.PlateCarree())
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES, edgecolor = 'blue')
#plt.plot(mesh_gridx, mesh_gridy, marker='o', color='k', linestyle='none')

#This is for adding radar data to the FED plot (optional)
#display.plot_ppi_map('reflectivity', sweep, vmin = -20, vmax=60, alpha = 0.1, cmap='pyart_HomeyerRainbow', colorbar_flag = False, title_flag = False, add_grid_lines = False, min_lat=radar_data.latitude['data'][0]-2, max_lat=radar_data.latitude['data'][0]+2, min_lon=radar_data.longitude['data'][0]-2, max_lon=radar_data.longitude['data'][0]+2, ax = ax)
lol = ax.imshow(save,vmin=0, vmax=60, cmap='pyart_HomeyerRainbow', alpha = 0.5, zorder = 10, extent=(grided_radar[0].origin_longitude['data'][0]-1.5, grided_radar[0].origin_longitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]-1.5))
ax.contour(fsave, vmin=0, vmax=np.max(fsave)/np.max(fsave), levels = levels, cmap='HomeyerRainbow', zorder = 15, alpha = 0.7, origin = 'upper', extent=(radar_frequency_files[0].origin_longitude['data'][0]-1.5, radar_frequency_files[0].origin_longitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]+1.5, radar_frequency_files[0].origin_latitude['data'][0]-1.5))
im = ax.pcolormesh(yedges, xedges, grid_points, norm = 'log', zorder = 11, cmap=load_cmap("LightBlue2DarkBlue10Steps"), vmin = 1, vmax = 1000, alpha = 1)
#ax.set_xlim(event_center[0]-zoom, event_center[0]+zoom)
#ax.set_ylim(event_center[1]-zoom, event_center[1]+zoom)
#ax.set_xlim(-77.160, -75.419)
ax.set_xlim(-77.7, -75.21)
#ax.set_ylim(43.33, 44.09)
ax.set_ylim(43.33, 44.2)
ax.plot(-76.506667, 43.454722, 'o', color = 'k', markersize = 7, zorder = 20)
ax.text(-76.496667, 43.464722, 'Oswego', color = 'k', fontsize = 12, zorder = 20)
ax.plot(-75.91142238037271, 43.97305891066422, 'o', color = 'k', markersize = 7, zorder = 20)
ax.text(-75.89142238037271, 43.99305891066422, 'Watertown', color = 'k', fontsize = 12, zorder = 20)
ax.plot(-75.68, 43.7558, 'o', color = 'k', markersize = 7, zorder = 20)
ax.text(-75.665, 43.7658, 'KTYX', color = 'k', fontsize = 12, zorder = 20)
ax.plot(-75.6364, 43.8218, 'o', color = 'k', markersize = 7, zorder = 20)
ax.text(-75.6464, 43.8178, 'Wind Turbine \n    Fields', color = 'k', fontsize = 12, zorder = 20)
#ax.axes('off')
#plt.imshow(grided_data, cmap = 'cool')
plt.tight_layout()
plt.show()

#%% n-Panel Plot of FED

image1 = Image.open('C:/Users/BenLa/Downloads/FED_0951-1214_PlotAvg_ContFrequ20_ver4.png')
image1_array = np.array(image1)

image2 = Image.open("C:/Users/BenLa/Downloads/FED_0815-1551_PlotAvg_ContFrequ20.png")
image2_array = np.array(image2)

image3 = Image.open("C:/Users/BenLa/Downloads/FED_2211-0443_PlotAvg_ContFrequ20_ver2.png")
image3_array = np.array(image3)

image4 = Image.open("C:/Users/BenLa/Downloads/FED_0738-1029_PlotAvg_ContFrequ20_ver4.png")
image4_array = np.array(image4)

image5 = Image.open('C:/Users/BenLa/Downloads/FED_1029-1402_plotAvg_ContFrequ20.png')
image5_array = np.array(image5)

image6 = Image.open('C:/Users/BenLa/Downloads/FED_1735-2117_PlotAvg_ContFrequ20_ver3.png')
image6_array = np.array(image6)

fig = plt.figure(figsize=(24, 15))

ax1 = fig.add_axes([0.01, 0.65, 0.39, 0.32])
ax1.axis('off')
fig.text(0.12, 0.93, '17 Nov 2022 0945-1215 UTC \n           30 Radar Scans', fontsize = 20)
fig.text(0.02, 0.89, '(a)', fontsize = 18, fontweight = 'bold')
ax2 = fig.add_axes([0.41, 0.65, 0.39, 0.32])
ax2.axis('off')
fig.text(0.52, 0.93, '18 Nov 2022 0815-1600 UTC \n           92 Radar Scans', fontsize = 20)
fig.text(0.42, 0.89, '(b)', fontsize = 18, fontweight = 'bold')
ax3 = fig.add_axes([0.01, 0.35, 0.39, 0.32])
ax3.axis('off')
fig.text(0.09, 0.63, '18 Nov 2022 2215 - 19 Nov 2022 0445 UTC \n                    72 Radar Scans', fontsize = 20)
fig.text(0.02, 0.59, '(c)', fontsize = 18, fontweight = 'bold')
ax4 = fig.add_axes([0.41, 0.35, 0.39, 0.32])
ax4.axis('off')
fig.text(0.52, 0.63, '20 Nov 2022 0730-1030 UTC \n           30 Radar Scans', fontsize = 20)
fig.text(0.42, 0.59, '(d)', fontsize = 18, fontweight = 'bold')
ax5 = fig.add_axes([0.01, 0.05, 0.39, 0.32])
ax5.axis('off')
fig.text(0.12, 0.33, '20 Nov 2022 1030-1402 UTC \n           36 Radar Scans', fontsize = 20)
fig.text(0.02, 0.29, '(e)', fontsize = 18, fontweight = 'bold')
ax6 = fig.add_axes([0.41, 0.05, 0.39, 0.32])
ax6.axis('off')
fig.text(0.52, 0.33, '20 Nov 2022 1735-2115 UTC \n           45 Radar Scans', fontsize = 20)
fig.text(0.42, 0.29, '(f)', fontsize = 18, fontweight = 'bold')

cax = fig.add_axes([0.82, 0.05, 0.03, 0.9])
cax3 = fig.add_axes([0.89, 0.05, 0.03, 0.9])
cax2 = fig.add_axes([0.89, 0.05, 0.03, 0.9])

ax1.imshow(image1_array)
ax2.imshow(image2_array)
ax3.imshow(image3_array)
ax4.imshow(image4_array)
ax5.imshow(image5_array)
ax6.imshow(image6_array)

cbar = fig.colorbar(im, cax, orientation='vertical')
cbar.set_label(label = '  LMA Points Per Grid Box', size = 20, rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=15)

levels2 = [6, 12, 18, 24, 30, 36, 42, 48, 54, 59.9]
cbar3 = fig.colorbar(ok, cax3, orientation='vertical', ticks = [0.19, 0.28, 0.37, 0.46, 0.55, 0.64, 0.73, 0.82, 0.91, 1])
cbar3.ax.set_yticklabels(['10%', '20%', '30%', '40%', '    ,  50%', '60%', '70%', '80%', '90%', '    ,  100%'])
cbar3.ax.tick_params(labelsize=15)
#cbar3.set_label(label = 'Frequency of Reflectivity > 20 dBZ per Grid Box', size = 20, rotation = 270, labelpad = 50)

cbar2 = fig.colorbar(lol, cax2, orientation='vertical')
cbar2.set_label(label = '  Mean Reflectivity in dBZ (Filled) \n Percent of Radar Scans with Reflectivity > 20 dBZ (Line)', size = 20, rotation=270, labelpad=95)
cbar2.add_lines(levels = levels2, colors = ['#026dc6', '#76bed1', '#7fcbb0', '#a0e185', '#ddfa67', '#e6da24', '#dfa700', '#d67000', '#ca3800', '#c42421'], linewidths = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
cbar2.ax.tick_params(labelsize=15)

plt.tight_layout()
plt.show()
