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
import pyxlma.plot.leader_speed as lmaleader
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


#%%

def negativeaxis(position, radar):
    
    if position < 0:
        if position < radar.longitude['data'][0]:
            return True
    else:
        if position < radar.latitude['data'][0]:
            return True
    return False

def VelocityPlot(interactive_lma):

        lat = interactive_lma.this_lma_lat
        lon = interactive_lma.this_lma_lon
        alt = interactive_lma.this_lma_alt
        time = interactive_lma.this_lma_time
        first = np.nanargmin(time)
        #print(time.min())
        #print(np.nanargmin(time))
        #print(time[first])
        

        distance_from_origin, time_from_origin = lmaleader.get_time_distance(lat, lon, time, lat[first], lon[first], time.min())
        
        fig = plt.figure(figsize = (10, 10))
        ax = fig.add_subplot()
        lmaleader.time_distance_plot(ax, time_from_origin, distance_from_origin, c = alt)
        
        plt.legend(title = "Leader Type Key", loc = 2, fontsize=12, title_fontsize=16, framealpha=1)
        
        plt.show()
                                                                                #mx, my, xx, xy
def RadPlanPlot(interactive_lma, radar, variable = 'reflectivity', xsec = [-76.26, 43.64, -75.43, 43.97, 15]):
        if type(radar) != pyart.core.Radar:
                radar =  pyart.io.read(radar)
        
        sweep=1 # sweep 0 = lowest level
        fig = plt.figure(figsize=(15,10))
        ax = plt.subplot(projection = ccrs.PlateCarree())
        display = pyart.graph.RadarMapDisplay(radar)
        display.plot_ppi_map(variable, sweep, vmin = -20, vmax=60, cmap='pyart_HomeyerRainbow', min_lat=radar.latitude['data'][0]-2, max_lat=radar.latitude['data'][0]+2, min_lon=radar.longitude['data'][0]-2, max_lon=radar.longitude['data'][0]+2, ax = ax)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
        ax.add_feature(cfeature.BORDERS)
        
        plt.plot([xsec[0], xsec[2]], [xsec[1], xsec[3]], color = '#6a8f9c', linewidth = 2)

        fruit = 0
        veggies = 0
        plop = []
        save_xsec_lma = []
        trash_xsec_lma = []
        max_dist_deg = 0.031
        start = interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[0]] - dt.timedelta(seconds = 2.5*60)
        end   = interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[-1]] + dt.timedelta(seconds = 7.5*60)

        while fruit < len(interactive_lma.this_lma_alt):
                if interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') >= start.strftime('%H:%M:%S.%f') and interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') <= end.strftime('%H:%M:%S.%f'):
                        #print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f'))
                        plop.append(np.cross((xsec[2], xsec[3])-np.array([xsec[0], xsec[1]]), np.array([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit]])-np.array([xsec[0], xsec[1]]))/norm((xsec[2], xsec[3])-np.array([xsec[0], xsec[1]])))
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
        
def RadXSecPlot(interactive_lma, radar, xsec = [-76.26, 43.64, -75.43, 43.97, 15]):
                
        if type(radar) != pyart.core.Radar:
                radar =  pyart.io.read(radar)
        
        radarx = radar.longitude['data'][0]
        radary = radar.latitude['data'][0]
        
        xsec0_rad, xsec2_rad = lmarad.haversine(xsec[1], xsec[0], xsec[1], radarx)/1e3, lmarad.haversine(xsec[3], xsec[2], xsec[3], radarx)/1e3
        xsec1_rad, xsec3_rad = lmarad.haversine(radary, xsec[0], xsec[1], xsec[0])/1e3, lmarad.haversine(radary, xsec[2], xsec[3], xsec[2])/1e3
        
        if negativeaxis(xsec[0], radar): xsec0_rad = xsec0_rad*-1
        if negativeaxis(xsec[1], radar): xsec1_rad = xsec1_rad*-1
        if negativeaxis(xsec[2], radar): xsec2_rad = xsec2_rad*-1
        if negativeaxis(xsec[3], radar): xsec3_rad = xsec3_rad*-1
        
        nn = int(((xsec2_rad-xsec0_rad)**2+(xsec3_rad-xsec1_rad)**2)**0.5/0.1) 
        des_x = np.linspace(xsec0_rad,xsec2_rad,nn)*1e3
        des_y = np.linspace(xsec1_rad,xsec3_rad,nn)*1e3
        des_z = np.arange(0,xsec[4]+0.1,0.1)*1e3

        desx_grid,desz_grid = np.meshgrid(des_x,des_z)
        desy_grid,desz_grid = np.meshgrid(des_y,des_z)
        new_x = np.arange(0,np.shape(des_x)[0]*0.1 ,0.1)*1e3

        keep_sweeps_pol  = [] # Holding spot for good polarimetric variable (non-NaN ZDR) sweeps
        keep_sweeps_trad = [] # Holding spot for good traditional variable (non-NaN velocity) sweeps
        for sweep in radar.sweep_number['data']:
            if not (np.size(radar.extract_sweeps([sweep]).fields['differential_reflectivity']['data'].mask)-
                     np.sum(radar.extract_sweeps([sweep]).fields['differential_reflectivity']['data'].mask))==0:
                #print ('Keeping polarimetric sweep number: ', sweep)
                keep_sweeps_pol+=[sweep]
            if not (np.size(radar.extract_sweeps([sweep]).fields['velocity']['data'].mask)-
                     np.sum(radar.extract_sweeps([sweep]).fields['velocity']['data'].mask))==0:
                #print ('Keeping traditional sweep number: ', sweep)
                keep_sweeps_trad+=[sweep]

        # Keep only the good sweeps for each set of radar variables
        radar_pol = radar.extract_sweeps(keep_sweeps_pol)
        radar_trad = radar.extract_sweeps(keep_sweeps_trad)

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
        
        fruit = 0
        veggies = 0
        plop = []
        save_xsec_lma = []
        trash_xsec_lma = []
        max_dist_deg = 0.031
        start = interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[0]] - dt.timedelta(seconds = 2.5*60)
        end   = interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[-1]] + dt.timedelta(seconds = 7.5*60)

        while fruit < len(interactive_lma.this_lma_alt):
                if interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') >= start.strftime('%H:%M:%S.%f') and interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f') <= end.strftime('%H:%M:%S.%f'):
                        #print(interactive_lma.this_lma_time[interactive_lma.this_lma_time.index[fruit]].strftime('%H:%M:%S.%f'))
                        plop.append(np.cross((xsec[2], xsec[3])-np.array([xsec[0], xsec[1]]), np.array([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit]])-np.array([xsec[0], xsec[1]]))/norm((xsec[2], xsec[3])-np.array([xsec[0], xsec[1]])))
                        #print(cake)
                        #print(d)
                        if abs(plop[veggies]) < max_dist_deg:
                                save_xsec_lma.append([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], interactive_lma.this_lma_alt[fruit]])
                        else:
                                trash_xsec_lma.append([interactive_lma.this_lma_lon[fruit], interactive_lma.this_lma_lat[fruit], interactive_lma.this_lma_alt[fruit]])
                        veggies += 1
                fruit += 1
        
        if 'kdpdata' not in locals():
                kdpdata = pyart.retrieve.kdp_vulpiani(radar_trad,  phidp_field ='differential_phase',  band ='S' , windsize = 34)
                kdpdata = kdpdata[0]['data'] #Retreive the kdp 'data' from the the kdpdata that was calculated using the Vulpiani method 
                mask = np.logical_and(kdpdata > -0.01, kdpdata < 0.01) #mask the values from -0.01 to 0.01 so they do not plot
                kdpdata = np.where(mask, np.nan, kdpdata) #Apply the mask to kdpdata
                radar_trad.add_field('specific_differential_phase_hv', {'data': kdpdata.data})  #Add a new field to the radar data dictionary with the kdp data
                
        lma_alt = interactive_lma.this_lma_alt
        new_lma_x = interactive_lma.this_lma_lon
        new_lma_y = interactive_lma.this_lma_lat

        t = 0

        new_lma_r = []
        new_lma_alt = []

        while t < len(save_xsec_lma)-1:
            new_lma_r.append(lmarad.haversine(xsec[0], xsec[1], save_xsec_lma[t][1], save_xsec_lma[t][0])/1e3)
            new_lma_alt.append(save_xsec_lma[t][2])
            #new_lma_r.append(haversine(my, mx, new_lma_y[t], new_lma_x[t])/1e3)
            t+=1
        fig = plt.figure(figsize=(28,16))
        #Reflectivity
        ax1 = fig.add_subplot(321)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['reflectivity']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(-20,72,1),cmap='pyart_HomeyerRainbow')
        plt.colorbar(label='Reflectivity (dBZ)')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        plt.ylim(0, xsec[4])
        #Velocity
        ax2 = fig.add_subplot(322)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['velocity']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(-40,40,1),cmap='NWSVel')
        plt.colorbar(label='Velocity (m/s)')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        plt.ylim(0, xsec[4])
        #Spectrum Width
        ax3 = fig.add_subplot(323)
        colors1 = plt.cm.binary_r(np.linspace(0.2,0.8,33))
        colors2 = plt.cm.gnuplot_r(np.linspace(0.,0.7,100))
        colors = np.vstack((colors1, colors2[10:121]))
        swcolours = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_trad.fields['spectrum_width']['data'].ravel()[indext_trad].reshape(np.shape(desz_grid)),levels = np.arange(0,14,0.1),cmap='NWS_SPW')
        plt.colorbar(label='Spectrum Width')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        #Differential Reflectivity
        ax4 = fig.add_subplot(324)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_pol.fields['differential_reflectivity']['data'].ravel()[indext_pol].reshape(np.shape(desz_grid)),levels = np.arange(-4,8,0.1),cmap='ChaseSpectral')
        plt.colorbar(label='Differential Reflectivity (ZDR)')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        plt.ylim(0, xsec[4])
        #Correlation Coefficient
        ax5 = fig.add_subplot(325)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, radar_pol.fields['cross_correlation_ratio']['data'].ravel()[indext_pol].reshape(np.shape(desz_grid)),levels = np.arange(0,1.1,0.01),cmap='SCook18')
        plt.colorbar(label='Correlation Coefficient')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        plt.ylim(0, xsec[4])
        #Specific Differential Phase
        ax6 = fig.add_subplot(326)
        plt.contourf(new_x/1e3, desz_grid[:,0]/1e3, kdpdata.ravel()[indext_pol].reshape(np.shape(desz_grid)),np.arange(-1,4.1,0.1),cmap=swcolours)
        plt.colorbar(label='Specific Differential Phase (KDP)')
        plt.scatter(new_lma_r, new_lma_alt, color='k', s=25)
        plt.xlim(0, np.max(new_x)/1e3)
        plt.ylim(0, xsec[4])
                 
                 
        plt.ylabel('Altitude (km AGL)')
        plt.xlabel('Distance (km)')
        plt.tight_layout()
        
def DOWPlot(interactive_lma, radar, dow, max_range = 50, max_z = 15):
        tbuffer = 2.5*60 # How many seconds before and after the radar scan time do you want to plot LMA data?
        max_dist_deg = 0.031 # degrees, max distance of an LMA source to a radar grid centroid to be plotted, roughtly 2.5km in midlatitudes

        # Set max range for DOW RHI
        rng_rhi = 55

        # Start with the DOW file of interest
        dow_file = dow
        # Read into pyart
        dow_pyart = pyart.io.read(dow_file)
        gatefilter = pyart.filters.GateFilter(dow_pyart)
        gatefilter.exclude_below('DBZHCC', 10)
        if type(radar) != pyart.core.Radar:
                radar =  pyart.io.read(radar)
        # Find x,y,z of gates
        rx,ry,rz = dow_pyart.get_gate_x_y_z(0)
        radarx = radar.longitude['data'][0]
        radary = radar.latitude['data'][0]
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
                                                            radar.longitude['data'][0], 
                                                            radar.latitude['data'][0], 
                                                            R=6370997.0)
        g2_x =  np.cos(az_rad)*rng_rhi
        g2_y = -np.sin(az_rad)*rng_rhi
        # And the lat/lon of that point
        END_lon, END_lat = pyart.core.cartesian_to_geographic_aeqd(g1_x+g2_x*1e3, g1_y+g2_y*1e3,
                                                                radar.longitude['data'][0], 
                                                                radar.latitude['data'][0], 
                                                                R=6370997.0)

        mx_DOW, xx_DOW = lmarad.haversine(DOW_lat, DOW_lon, DOW_lat, radarx)/1e3, lmarad.haversine(END_lat, END_lon, END_lat, radarx)/1e3
        my_DOW, xy_DOW = lmarad.haversine(radary, DOW_lon, DOW_lat, DOW_lon)/1e3, lmarad.haversine(radary, END_lon, END_lat, END_lon)/1e3

        nn = int(((xx_DOW-mx_DOW)**2+(xy_DOW-my_DOW)**2)**0.5/0.1) 
        DOW_des_x = np.linspace(mx_DOW,xx_DOW,nn)*1e3
        DOW_des_y = np.linspace(my_DOW,xy_DOW,nn)*1e3
        DOW_des_z = np.arange(0,max_z+0.1,0.1)*1e3

        desx_grid,desz_grid = np.meshgrid(DOW_des_x,DOW_des_z)
        desy_grid,desz_grid = np.meshgrid(DOW_des_y,DOW_des_z)
        DOWnew_x = np.arange(0,np.shape(DOW_des_x)[0]*0.1 ,0.1)*1e3

        sweep=1 # sweep 0 = lowest level
        fig = plt.figure(figsize=(15,10))
        ax = plt.subplot(projection = ccrs.PlateCarree())
        display = pyart.graph.RadarMapDisplay(radar)
        display.plot_ppi_map('reflectivity', sweep, vmin = -20, vmax=60, cmap='pyart_HomeyerRainbow', min_lat=radar.latitude['data'][0]-2, max_lat=radar.latitude['data'][0]+2, min_lon=radar.longitude['data'][0]-2, max_lon=radar.longitude['data'][0]+2, ax = ax)
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
        
        wheat = 0
        DOW_new_lma_r = []
        DOW_lma_alt = []
        while wheat < len(save_lma):
            DOW_new_lma_r.append(lmarad.haversine(DOW_lat, DOW_lon, save_lma[wheat][1], save_lma[wheat][0])/1e3)
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

def AvgRad(radlist):
        radar_pickle = radlist[0]+str(len(radlist))+'.pickle'
        reflectivity_grids = []

        if os.path.exists(radar_pickle):
            with open(radar_pickle, 'rb') as pickle_file:
                grided_radar = pickle.load(pickle_file) 
        else:
                grided_radar, start_times = lmarad.ReadRadar(radlist)

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
        return save

def FrequRad(radlist):
        radar_pickle = radlist[0]+str(len(radlist))+'.pickle'
        levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        frequency_grids = []

        if os.path.exists(radar_pickle):
            with open(radar_pickle, 'rb') as pickle_file:
                radar_frequency_files = pickle.load(pickle_file) 
        else:
                radar_frequency_files, radar_frequency_start = lmarad.ReadRadar(radlist)

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
        
        return fsave

def FlashExtentDensity(interactive_lma, radlist, avg, frequ):
        
        radar_pickle = radlist[0]+str(len(radlist))+'.pickle'
        
        levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        if os.path.exists(radar_pickle):
            with open(radar_pickle, 'rb') as pickle_file:
                grided_radar = pickle.load(pickle_file) 
        else:
                grided_radar, start_times = lmarad.ReadRadar(radlist)

                with open(radar_pickle, 'wb') as pickle_file:
                        pickle.dump(grided_radar, pickle_file)
                        
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
        ax = grid_fig.add_axes([0.02, 0.05, 0.78, 0.9], projection = ccrs.PlateCarree())
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.STATES, edgecolor = 'blue')
        #plt.plot(mesh_gridx, mesh_gridy, marker='o', color='k', linestyle='none')

        #This is for adding radar data to the FED plot (optional)
        #display.plot_ppi_map('reflectivity', sweep, vmin = -20, vmax=60, alpha = 0.1, cmap='pyart_HomeyerRainbow', colorbar_flag = False, title_flag = False, add_grid_lines = False, min_lat=radar_data.latitude['data'][0]-2, max_lat=radar_data.latitude['data'][0]+2, min_lon=radar_data.longitude['data'][0]-2, max_lon=radar_data.longitude['data'][0]+2, ax = ax)
        lol = ax.imshow(avg,vmin=0, vmax=60, cmap='pyart_HomeyerRainbow', alpha = 0.5, zorder = 10, extent=(grided_radar[0].origin_longitude['data'][0]-1.5, grided_radar[0].origin_longitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]-1.5))
        ok = ax.contour(frequ, vmin=0, vmax=np.max(frequ)/np.max(frequ), levels = levels, cmap='HomeyerRainbow', zorder = 15, alpha = 0.7, origin = 'upper', extent=(grided_radar[0].origin_longitude['data'][0]-1.5, grided_radar[0].origin_longitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]+1.5, grided_radar[0].origin_latitude['data'][0]-1.5))
        im = ax.pcolormesh(yedges, xedges, grid_points, norm = 'log', zorder = 11, cmap=load_cmap("LightBlue2DarkBlue10Steps"), vmin = 1, vmax = 1000, alpha = 1)
        #ax.set_xlim(event_center[0]-zoom, event_center[0]+zoom)
        #ax.set_ylim(event_center[1]-zoom, event_center[1]+zoom)
        
        cax = grid_fig.add_axes([0.82, 0.05, 0.03, 0.9])
        cax3 = grid_fig.add_axes([0.89, 0.05, 0.03, 0.9])
        cax2 = grid_fig.add_axes([0.89, 0.05, 0.03, 0.9])
        
        cbar = grid_fig.colorbar(im, cax, orientation='vertical')
        cbar.set_label(label = '', size = 20, rotation=270, labelpad=15)
        cbar.ax.tick_params(labelsize=10)
        levels2 = [6, 12, 18, 24, 30, 36, 42, 48, 54, 59.9]
        cbar3 = grid_fig.colorbar(ok, cax3, orientation='vertical', ticks = [0.19, 0.28, 0.37, 0.46, 0.55, 0.64, 0.73, 0.82, 0.91, 1])
        cbar3.ax.set_yticklabels(['10%', '20%', '30%', '40%', '    ,  50%', '60%', '70%', '80%', '90%', '    ,  100%'])
        cbar3.ax.tick_params(labelsize=10)
        cbar2 = grid_fig.colorbar(lol, cax2, orientation='vertical')
        cbar2.set_label(label = '', size = 20, rotation=270, labelpad=95)
        cbar2.add_lines(levels = levels2, colors = ['#026dc6', '#76bed1', '#7fcbb0', '#a0e185', '#ddfa67', '#e6da24', '#dfa700', '#d67000', '#ca3800', '#c42421'], linewidths = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
        cbar2.ax.tick_params(labelsize=10)
        
        #ax.set_xlim(-77.160, -75.419)
        ax.set_xlim(-77.33, -75.21)
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











