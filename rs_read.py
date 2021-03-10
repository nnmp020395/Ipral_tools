"""
1st reading of two radiosondes files:

/bdd/GRUAN/RADIOSONDAGE/L2/TRP/2020/2020_12_09/TRP-RS-01_2_M10-GDP-BETA_001_20201209T120000_1-000-001.nc
/bdd/SIRTA/pub/basesirta/1a/rs/2019/08/29/rs_1a_LtrappesPptuv_v01_20190829_231600_120.nc
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

from pathlib import Path
import sys
import xarray as xr

import numpy as np
import pandas as pd
import datetime as dt

def read_rs_SIRTA(rs_file_path, outdir):
    """
    #__ before 31.08.2019 -> output = dataframe with columns: ['alt', 'lat', 'lon', 'ta', 'pres', 'rh']
    Input: 
    --rs_file_path: full path of rs file
    --outdir: path of directory for output image
    
    Output: 
    --plot which draws Pression, Temperature, RH and Altitude above sea level
    --simple netCDF 
    """
    # rs_file_path = Path('/bdd/SIRTA/pub/basesirta/1a/rs/2019/08/29/rs_1a_LtrappesPptuv_v01_20190829_231600_120.nc')
    rs_data = xr.open_dataset(rs_file_path)
    df = rs_data.to_dataframe()[['alt', 'lat', 'lon', 'ta', 'pres', 'rh']]
    resolution = np.float(rs_data.attrs['geospatial_vertical_resolution'].split(' ')[0])
    #__ 
    f, ax = plt.subplots(1,3, figsize=(15,10))
    ax[0].plot(df['pres'], df['alt'])
    ax[0].set(xlabel = rs_data.pres.standard_name + '_' + rs_data.pres.units, 
        ylabel='altitude_' + rs_data.alt.units +'\nvertical_resolution='+rs_data.attrs['geospatial_vertical_resolution'])
    ax[1].plot(df['ta'], df['alt'])
    ax[1].set(xlabel = rs_data.ta.standard_name + '_' + rs_data.ta.units)
    ax[2].plot(df['rh'], df['alt'])
    ax[2].set(xlabel = rs_data.rh.standard_name + '_' + rs_data.rh.units)
    plt.suptitle('Radiosondes MeteoFrance - SIRTA \n time_coverage_start '+rs_data.attrs['time_coverage_start'], fontweight='bold')
    plt.savefig(Path(outdir, rs_file_path.name.split('.')[0]+'PTRh.png'))


def read_rs_GRUAN(rs_file_path, outdir):
    """
    #__ after 31.08.2019 

    Input: 
    --rs_file_path: full path of rs file
    --outdir: path of directory for output image
    
    Output: 
    --plot which draws Pression, Temperature, RH and Altitude above sea level
    --simple netCDF 
    """
    # RS_PATH = Path('/bdd/GRUAN/RADIOSONDAGE/L2/TRP/2020/2020_12_09/TRP-RS-01_2_M10-GDP-BETA_001_20201209T120000_1-000-001.nc')
    rs_data = xr.open_dataset(rs_file_path)
# keys = list(rs_data.keys())
# for k in keys:
#     if (k.split('_')[0] == 'rh'): #|(k.split('_')[0] == 'press'):
#         print(f'variables: {k}')
#         print(rs_data[k].long_name)
#         print(rs_data[k].coords)
    time = rs_data.time.values
    alt_raw = rs_data.alt_raw.to_dataframe()
    press = rs_data.press.to_dataframe()
    temp = rs_data.temp.to_dataframe()
    rh_raw = rs_data.rh_raw.to_dataframe()
    lat = rs_data.lat.to_dataframe()
    lon = rs_data.lon.to_dataframe()
    resolution = np.float(rs_data.attrs['geospatial_vertical_resolution'].split(' ')[0])
    df = alt_raw.join([lat, lon, temp, press, rh_raw], how='left')
    #__ 
    f, ax = plt.subplots(1,3, figsize=(15,10))
    ax[0].plot(df['press'], df['alt_raw'])
    ax[0].set(xlabel = rs_data.press.standard_name + '_' + rs_data.press.units, 
        ylabel='altitute_' + rs_data.alt_raw.units +'\nvertical_resolution='+rs_data.attrs['geospatial_vertical_resolution'])
    ax[1].plot(df['temp'], df['alt_raw'])
    ax[1].set(xlabel = rs_data.temp.standard_name + '_' + rs_data.temp.units)
    ax[2].plot(df['rh_raw'], df['alt_raw'])
    ax[2].set(xlabel = rs_data.rh_raw.standard_name + '_' + rs_data.rh_raw.units)
    plt.suptitle('Radiosondes MeteoFrance - SIRTA \n time_coverage_start '+rs_data.attrs['time_coverage_start'], fontweight='bold')
    plt.savefig(Path(outdir, rs_file_path.name.split('.')[0]+'PTRh.png'))


from argparse import Namespace, ArgumentParser
parser = ArgumentParser()
parser.add_argument('--year', '-y', type=str, required=True)
parser.add_argument('--month', '-m',type=str, required=True)
parser.add_argument('--day', '-d',type=str, required=True)
parser.add_argument('--outdir', '-o',type=str, required=True)
opts = parser.parse_args()
#__list rs files 
year_str = opts.year
month_str = opts.month
day_str = opts.day
outdir = opts.outdir #Path('/homedata/nmpnguyen/IPRAL/1a_Fig')
if dt.date(int(year_str),int( month_str), int(day_str)) < dt.date(2019, 8, 31):
    MAIN_PATH = Path('/bdd/SIRTA/pub/basesirta/1a/rs/')
    RS_PATH = Path(MAIN_PATH, year_str, month_str, day_str)
    RS_MASK = 'rs_1a_LtrappesPptuv_v01_*.nc'
    rs_files = sorted(RS_PATH.glob(RS_MASK))
else:
    MAIN_PATH = Path('/bdd/GRUAN/RADIOSONDAGE/L2/TRP/')
    date_str = '_'.join([year_str, month_str, day_str])
    RS_PATH = Path(MAIN_PATH, year_str, date_str)
    RS_MASK = 'TRP-RS-01_2_M10-GDP-*.nc'
    rs_files = sorted(RS_PATH.glob(RS_MASK))


if len(rs_files) != 0:
    print(f'List of files found: \n{rs_files}')
else:
    print('No file found')


if dt.date(int(year_str),int( month_str), int(day_str)) < dt.date(2019, 8, 31): 
    for rs_file_path in rs_files:
        read_rs_SIRTA(rs_file_path, outdir)   
else:
    for rs_file_path in rs_files:
        read_rs_GRUAN(rs_file_path, outdir)


#__function get Pression or Temperature from radiosundes data
