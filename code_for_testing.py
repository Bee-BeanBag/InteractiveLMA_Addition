#%% Imports

import pyart
import urllib.request
import zipfile
import xarray as xr
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
import InteractiveLMA_Addition.Thesis_Code_Functs as lmafuncts
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
#%% LMA Setup

obsticles = "https://aeronav.faa.gov/Obst_Data/DOF_250511.zip"
urllib.request.urlretrieve(obsticles, "DOF_250511.zip")

zipdata = zipfile.ZipFile('DOF_250511.zip')
zipinfos = zipdata.infolist()

for zipinfo in zipinfos:
    # This will do the renaming
    if '.csv' in zipinfo.filename:
        print("Found ", zipinfo.filename)
    # We can override the filename. Neat!\n",
    zipinfo.filename = 'turbine_locations.csv'
    zipdata.extract(zipinfo)
    
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
        
ds = xr.open_dataset('./LYLOUT_221117_080000_360000_map4000NovFull1.nc')
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
DOF = './NewYork_Obsticles.Dat'

specs = [(35, 37), (37, 41), (41, 46), (48, 52), (52, 55), (55, 60), (82, 88), (62, 74)]
StateObs2 = pd.read_fwf(DOF, colspecs = specs, skiprows=4, names=('ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec', 'Elevation (ft)', 'Obsticle'))
#NYObs2 = pd.read_fwf('D:/NewYork_Obsticles.Dat', colspecs = specs, skiprows=9209, skipfooter =1,  names=('ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec', 'Elevation (ft)', 'Obsticle'))

#Bee's
#convert the deg/min/sec for lat and long to decimal format
StateObs2['ylat']=StateObs2['ylatdeg']+StateObs2['ylatmin']/60+StateObs2['ylatsec']/3600
StateObs2['xlong']=(StateObs2['xlongdeg']+StateObs2['xlongmin']/60+StateObs2['xlongsec']/3600)*-1

#Bee's
#get rid of the deg/min/sec columns and send it to a csv file
StateObs2.drop(columns=['ylatdeg', 'ylatmin', 'ylatsec', 'xlongdeg', 'xlongmin', 'xlongsec']).to_csv(r'./State_Obsticles.csv')
StateObs3 = pd.read_csv('./State_Obsticles.csv')

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
        

xsecLatLonAlt = [-76.26, 43.64, -75.43, 43.97, 15]

tlim = pd.to_datetime(starttime).to_pydatetime()-pd.Timedelta(hours=1), pd.to_datetime(endtime).to_pydatetime()-pd.Timedelta(hours=1)

#%% LMA Plot
interactive_lma = AnnotatedLMAPlot(ds, clon=lma_ctr_lon, clat=lma_ctr_lat, tlim=tlim) #, network_data = entln_data)

#interactive_lma = SuperLMAPlot(ds, clon=lma_ctr_lon, clat=lma_ctr_lat,tlim=tlim, network_data = entln_data, radar_data = radar_data, points = xsecLatLonAlt)

#%%
#lmafuncts.VelocityPlot(interactive_lma)

#lmafuncts.RadPlanPlot(interactive_lma, "./KTYX20221120_115435_V06", xsec = [-77.11, 43.41, -75.68, 43.7, 15])

#lmafuncts.RadXSecPlot(interactive_lma, "./KTYX20221120_115435_V06", xsec = [-77.11, 43.41, -75.68, 43.7, 15])

#lmafuncts.DOWPlot(interactive_lma, "./KTYX20221120_115435_V06", "./cfrad.20221120_115541.398_DOW7_v155_s00_az9.90_RHI")

radar_list = ["./KTYX20221120_115435_V06"]
avg = lmafuncts.AvgRad(radar_list)
frequ = lmafuncts.FrequRad(radar_list)

lmafuncts.FlashExtentDensity(interactive_lma, radar_list, avg, frequ)

