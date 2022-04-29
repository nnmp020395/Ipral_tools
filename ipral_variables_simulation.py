import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import xarray as xr
from pathlib import Path
import sys
# from argparse import Namespace, ArgumentParser
# import click
import numpy as np
import pandas as pd
from tqdm import tqdm
__version__ = "0.1.0"


# parser = ArgumentParser()
# parser.add_argument("--ipral_path", "-ip", type = str, required=True)
# opts = parser.parse_args()
# # path_raw = Path("/homedata/nmpnguyen/IPRAL/1a/2019-05-14_ipral_1a_bck_corrected.nc")
# ipral_file = Path(opts.ipral_path)


def variables_from_era(ipral_file):
    """
    Le script permet de lire input ERA5 et outpt Pression/Temperature
    input = ipral file path
    
    """
    print('-----GET IPRAL BCK CORRECTED FILE-----')
    d = xr.open_dataset(ipral_file)
    time = d.time.values
    YEAR = pd.to_datetime(time[0]).strftime('%Y')
    MONTH = pd.to_datetime(time[0]).strftime('%m')
    lon_ipral = round(4*float(d.geospatial_lon_min))/4 #round(float(d.geospatial_lon_min),2)
    lat_ipral = round(4*float(d.geospatial_lat_min))/4
    print(f'longitude: {lon_ipral}')
    print(f'latitude: {lat_ipral}')
    #----
    print('-----GET ERA5 FILE-----')
    ERA_FOLDER = Path("/bdd/ERA5/NETCDF/GLOBAL_025/hourly/AN_PL")

    ERA_FILENAME = YEAR+MONTH+".ap1e5.GLOBAL_025.nc"
    GEOPT_PATH = ERA_FOLDER / YEAR / Path("geopt."+ERA_FILENAME)
    TA_PATH = ERA_FOLDER / YEAR / Path("ta."+ERA_FILENAME)
    geopt = xr.open_dataset(GEOPT_PATH)
    ta = xr.open_dataset(TA_PATH)
    #----
    print('-----CONVERT TIME AND LOCALISATION-----')
    # date_start = pd.to_datetime(time[0])
    # date_end = pd.to_datetime(time[-1])
    time = pd.to_datetime(time).strftime('%Y-%m-%dT%H:00:00.000000000')
    time = time.astype('datetime64[ns]')
    time_unique = np.unique(time)
    LAT = geopt.latitude[np.where(np.abs(geopt.latitude.values - lat_ipral) <=0.25)[0][1]].values
    LON = geopt.longitude[np.where(np.abs(geopt.longitude.values - lon_ipral) <=0.25)[0][1]].values
    #----
    if pd.to_datetime(time[-1]).strftime('%m')==pd.to_datetime(time[0]).strftime('%m'):
        print(f'path of temperature {TA_PATH}')
        print(f'path of geopotential {GEOPT_PATH}')
    else:
        YEAR = pd.to_datetime(time[-1]).strftime('%Y')
        MONTH = pd.to_datetime(time[-1]).strftime('%m')
        ERA_FILENAME = YEAR+MONTH+".ap1e5.GLOBAL_025.nc"
        GEOPT_PATH = [GEOPT_PATH, ERA_FOLDER / YEAR / Path("geopt."+ERA_FILENAME)]
        TA_PATH = [TA_PATH, ERA_FOLDER / YEAR / Path("ta."+ERA_FILENAME)]
        print(f'path of temperature {TA_PATH}')
        print(f'path of geopotential {GEOPT_PATH}')
        geopt1 = xr.open_dataset(GEOPT_PATH[1])
        ta1 = xr.open_dataset(TA_PATH[1])
        geopt = xr.concat([geopt, geopt1], dim='time')
        ta = xr.concat([ta, ta1], dim='time')

    print(geopt, ta)
    from timeit import default_timer as timer
    TIME = timer()    
    geopt_for_ipral = geopt.sel(time=time_unique, latitude=LAT, longitude=LON)#.to_dataframe()#['geopt']
    ta_for_ipral = ta.sel(time=time_unique, latitude=LAT, longitude=LON)#.to_dataframe()#['ta']
    print(geopt_for_ipral, ta_for_ipral)
    print(f'Time loading {timer()-TIME}')
    #----
    print('-----GETTING PRESSURE AND TEMPERATURE-----')
    lat_ipral = np.deg2rad(lat_ipral)
    acc_gravity = 9.78032*(1+5.2885e-3*(np.sin(lat_ipral))**2 - 5.9e-6*(np.sin(2*lat_ipral))**2)
    r0 = 2*acc_gravity/(3.085462e-6 + 2.27e-9*np.cos(2*lat_ipral) - 2e-12*np.cos(4*lat_ipral))
    g0 = 9.80665
    geopt_for_ipral['geopt_height'] = geopt_for_ipral["geopt"]/g0
    geopt_for_ipral['altitude'] = (geopt_for_ipral['geopt_height']*r0)/(acc_gravity*r0/g0 - geopt_for_ipral['geopt_height'])
    M = 28.966E-3 
    R = 8.314510
    T = (15 + 273.15)
    const = -(M*g0)/(R*T)
    p0 = 101325
    geopt_for_ipral['pression'] = p0*np.exp(const*geopt_for_ipral['altitude'])
    output_era = pd.merge(geopt_for_ipral, ta_for_ipral['ta'], left_index=True, right_index=True) 
    # output_era = xr.merge([geopt_for_ipral, ta_for_ipral])
    # output_era = output_era.to_dataframe()
    print('variables_from_era --> end')
    return output_era

def simulate_atb_mol(era):
    """
    Input is the output dataframe of variables_from_era function
    """
    print('-----SIMULATE ATTENUATED BACKSCATTERING FROM ERA5-----')
    k = 1.38e-23
    const355 = (5.45e-32/1.38e-23)*((355e-3/0.55)**(-4.09))
    const532 = (5.45e-32/1.38e-23)*((532e-3/0.55)**(-4.09))
    era['beta355'] = const355*era['pression'].div(era['ta'])
    era['beta532'] = const532*era['pression'].div(era['ta'])
    ratio = 0.119/3
    era['alpha355'] = era['beta355']/ratio
    era['alpha532'] = era['beta532']/ratio
    era = era.sort_index()
    level = np.unique(era.index.get_level_values(0))
    time = np.unique(era.index.get_level_values(1)) 
    era['tau355']=0 ;  era['tau532'] = 0
    A = pd.DataFrame()
    for t in time:
        a = era.loc[pd.IndexSlice[:,t],:].sort_index(ascending = False)
        for i in range(1, a.shape[0]):
            a['tau355'].iloc[i] = a['tau355'].iloc[i-1] + a['alpha355'].iloc[i]*(a['altitude'].iloc[i]-a['altitude'].iloc[i-1])
            a['tau532'].iloc[i] = a['tau532'].iloc[i-1] + a['alpha532'].iloc[i]*(a['altitude'].iloc[i]-a['altitude'].iloc[i-1])
        A = pd.concat((A, a), axis=0)
    A['beta355mol'] = A['beta355']*np.exp(-2*A['tau355'])
    A['beta532mol'] = A['beta532']*np.exp(-2*A['tau532'])
    print('simulate_atb_mol --> end')
    return A


import scipy.interpolate as spi
def interpolate_atb_mol(ipral_file, era, output):     
    """
    the Input is the output dataframe of simulate_atb_mol function
    """
    print('-----BEFORE INTERPOLATE-----')
    d = xr.open_dataset(ipral_file)
    r = d.range.values
    timeIpral = d.time.values
    timeEra = np.unique(era.index.get_level_values(1)) 
    time_tmp = np.array(pd.to_datetime(timeIpral).strftime('%Y-%m-%dT%H:00:00')).astype('datetime64[ns]')
    if len(time_tmp) != len(timeIpral):
        print("Time Error")
        sys.exit(1)
    #------
    columns_names = ['altitude', 'pression', 'ta', 'beta355mol', 'beta532mol']#, 'beta355', 'beta532', 'alpha355', 'alpha532', 'tau355', 'tau532'
    beta355mol_interp ,beta532mol_interp, pression_interp, ta_interp = [[] for _ in range(len(columns_names)-1)] #beta355_interp ,beta532_interp ,tau355_interp ,tau532_interp ,alpha355_interp ,alpha532_interp ,
    new_index = pd.MultiIndex.from_product([timeIpral, r], names = ['time', 'range'])
    # df_new = pd.DataFrame(index = new_index, columns = era.columns)
    print('-----INTERPOLATE ATTENUATED BACKSCATTERING FROM ERA5-----')
    for t1 in tqdm(time_tmp):
        a = era.loc[pd.IndexSlice[:, t1], columns_names]
        f1 = spi.interp1d(a['altitude'], a['beta355mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f2 = spi.interp1d(a['altitude'], a['beta532mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f3 = spi.interp1d(a['altitude'], a['tau355'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f4 = spi.interp1d(a['altitude'], a['tau532'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f5 = spi.interp1d(a['altitude'], a['alpha355'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f6 = spi.interp1d(a['altitude'], a['alpha532'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f7 = spi.interp1d(a['altitude'], a['beta355'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f8 = spi.interp1d(a['altitude'], a['beta532'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f9 = spi.interp1d(a['altitude'], a['pression'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f10 = spi.interp1d(a['altitude'], a['ta'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        beta355mol_interp, beta532mol_interp = np.append(beta355mol_interp, np.array(f1(r))), np.append(beta532mol_interp, np.array(f2(r)))
        # tau355_interp, tau532_interp = np.append(tau355_interp, np.array(f3(r))), np.append(tau532_interp, np.array(f4(r)))
        # alpha355_interp, alpha532_interp = np.append(alpha355_interp, np.array(f5(r))), np.append(alpha532_interp, np.array(f6(r)))
        # beta355_interp, beta532_interp = np.append(beta355_interp, np.array(f7(r))), np.append(beta532_interp, np.array(f8(r)))
        pression_interp, ta_interp = np.append(pression_interp, np.array(f9(r))), np.append(ta_interp, np.array(f10(r)))
        # print(str(t1))
    #------
    new_df = pd.DataFrame(index = new_index, data = np.array([pression_interp, ta_interp ,beta355mol_interp ,beta532mol_interp]).T, columns = columns_names[1:])
    #, beta355_interp ,beta532_interp ,alpha355_interp ,alpha532_interp ,tau355_interp ,tau532_interp
    #new_df.to_pickle("/homedata/nmpnguyen/IPRAL/1a/"+ipral_file.name.split(".")[0]+"_simul.pkl")
    new_df.to_pickle(Path(output, ipral_file.name.split(".")[0]+"_simul.pkl"))
    print('interpolate_atb_mol --> end')
    print(f"output filename: {Path(output, ipral_file.name.split('.')[0]+'_simul.pkl')}")
    
    return new_df


def simulate(ipral_file):
    # if option == "era": 
    out = variables_from_era(ipral_file)
    era = simulate_atb_mol(out)
    new_df = interpolate_atb_mol(ipral_file, era, output=Path("/homedata/nmpnguyen/IPRAL/RF/Simul/"))
    return new_df
    # else:
    #     return 0


#ipral_file = sorted(Path('/bdd/SIRTA/pub/basesirta/1a/ipral/2020/02/06/').glob('**/*_Lz1R15mF30sPbck_v01_*.nc'))
### Ouvrir le fichier des jours sélectionés
# with open('/homedata/nmpnguyen/ipral-tools/ClearSkyIprallist.txt', 'r') as f:
#     all_data = [line.strip() for line in f.readlines()]
    
# metadata_line = all_data[:4]
# listdays = all_data[4:]

# for l in listdays:
#   pf = list(Path('/bdd/SIRTA/pub/basesirta/1a/ipral/', l).glob('ipral_1a_Lz1R15mF30sPbck_v01_*.nc'))[0]
#   simulate(pf, 'era')

#__ get rs path corresponding to ipral data __
"""
ipral path --> get day --> if date before 31/08/2019 get path from baseipral/1a/rs, if after get path from GRUAN directory
sub-function: read data from GRUAN or BASESIRTA --> output = dataframe which store 'lat', 'lon', 'alt', 'ta', 'pression'

ignore time 
"""
