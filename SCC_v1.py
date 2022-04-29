import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

from pathlib import Path
import sys
import xarray as xr

import numpy as np
import pandas as pd

sys.path.append('/homedata/nmpnguyen/ipral-tools/')
from ipral_variables_simulation import simulate_atb_mol
"""
Objectif:   ouvrir les fichiers SCC produits de Chrstophes pour calibrer 
            comparer avec les ipral qui sont calibrées de meme facon

"""

list_scc = sorted(Path('/homedata/pietras/IPRAL/SCC_Products/2020/2020_09_14/Near_Range/').glob('20200914*.nc'))

for l in range(len(list_scc)):
    fl = list_scc[l]
    scc = xr.open_dataset(fl)
    alt = scc['altitude'].values[0,:]
    time = scc['time'].values
    YEAR = pd.to_datetime(time[0]).strftime('%Y')
    MONTH = pd.to_datetime(time[0]).strftime('%m')
    lon_ipral = float(scc['longitude'])#round(float(d.geospatial_lon_min),2)
    lat_ipral = float(scc['latitude'])

    print('-----GET ERA5 FILE-----')
    ERA_FOLDER = Path("/bdd/ERA5/NETCDF/GLOBAL_025/hourly/AN_PL")
    ERA_FILENAME = YEAR+MONTH+".ap1e5.GLOBAL_025.nc"
    GEOPT_PATH = ERA_FOLDER / YEAR / Path("geopt."+ERA_FILENAME)
    TA_PATH = ERA_FOLDER / YEAR / Path("ta."+ERA_FILENAME)
    print(f'path of temperature {TA_PATH}')
    print(f'path of geopotential {GEOPT_PATH}')
    geopt = xr.open_dataset(GEOPT_PATH)
    ta = xr.open_dataset(TA_PATH)
    print('-----CONVERT TIME-----')
    # date_start = pd.to_datetime(time[0])
    # date_end = pd.to_datetime(time[-1])
    time = pd.to_datetime(time).strftime('%Y-%m-%dT%H:00:00.000000000')
    time = time.astype('datetime64[ns]')
    time_unique = np.unique(time)

    LAT = geopt.latitude[np.where(np.abs(geopt.latitude.values - lat_ipral) <=0.25)[0][1]].values
    LON = geopt.longitude[np.where(np.abs(geopt.longitude.values - lon_ipral) <=0.25)[0][1]].values
    geopt_for_ipral = geopt.sel(time=time_unique, latitude=LAT, longitude=LON).to_dataframe()#['geopt']
    ta_for_ipral = ta.sel(time=time_unique, latitude=LAT, longitude=LON).to_dataframe()#['ta']
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
    print('variables_from_era --> end')
    output_era_allparam = simulate_atb_mol(output_era)

    timeSCC = scc['time'].values
    timeEra = np.unique(output_era_allparam.index.get_level_values(1)) 
    time_tmp = pd.to_datetime(timeSCC).strftime('%Y-%m-%dT%H:00:00.000000000').astype('datetime64[ns]')
    if len(time_tmp) != len(timeSCC):
        print("Time Error")
        sys.exit(1)
        #------
    columns_names = ['altitude', 'pression', 'ta', 'beta355mol', 'beta532mol']#, 'beta355', 'beta532', 'alpha355', 'alpha532', 'tau355', 'tau532'
    beta355mol_interp ,beta532mol_interp, pression_interp, ta_interp = [[] for _ in range(len(columns_names)-1)] #beta355_interp ,beta532_interp ,tau355_interp ,tau532_interp ,alpha355_interp ,alpha532_interp ,
    new_index = pd.MultiIndex.from_product([timeSCC, alt], names = ['time', 'alt'])
    # df_new = pd.DataFrame(index = new_index, columns = era.columns)
    print('-----INTERPOLATE ATTENUATED BACKSCATTERING FROM ERA5-----')

    import scipy.interpolate as spi    
    for t1 in time_tmp:
        a = output_era_allparam.loc[pd.IndexSlice[:, t1], columns_names]
        f1 = spi.interp1d(a['altitude'], a['beta355mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f2 = spi.interp1d(a['altitude'], a['beta532mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f9 = spi.interp1d(a['altitude'], a['pression'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f10 = spi.interp1d(a['altitude'], a['ta'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        beta355mol_interp, beta532mol_interp = np.append(beta355mol_interp, np.array(f1(alt))), np.append(beta532mol_interp, np.array(f2(alt)))
        pression_interp, ta_interp = np.append(pression_interp, np.array(f9(alt))), np.append(ta_interp, np.array(f10(alt)))
    #------
    new_df = pd.DataFrame(index = new_index, data = np.array([pression_interp, ta_interp ,beta355mol_interp ,beta532mol_interp]).T, columns = columns_names[1:])
    # new_df.to_pickle("/homedata/nmpnguyen/IPRAL/1a/"+ipral_file.name.split(".")[0]+"_simul.pkl")
    print('interpolate_atb_mol --> end')

    beta532mol = new_df['beta532mol']
    beta532mol = beta532mol.unstack(level=1).values
    beta355mol = new_df['beta355mol']
    beta355mol = beta355mol.unstack(level=1).values
    print('SIMULATED DATA-------------------------------end')

    rcs0 = scc['range_corrected_signal'][0,:,:].values
    rcs1 = scc['range_corrected_signal'][1,:,:].values
    atm0 = scc['atmospheric_background'][0,:].values
    atm1 = scc['atmospheric_background'][1,:].values
    scc0 = rcs0 #- atm0.reshape((atm0.size,1))
    scc1 = rcs1 #- atm1.reshape((atm1.size,1))
    scc0Av = pd.DataFrame(scc0).groupby(np.arange(scc0.shape[1])//8, axis=1).mean().values
    scc1Av = pd.DataFrame(scc1).groupby(np.arange(scc1.shape[1])//8, axis=1).mean().values
    z_cc = np.where((alt[::8]>=4000)&(alt[::8]<=4200))[0]
    beta355molAv = pd.DataFrame(beta355mol).groupby(np.arange(beta355mol.shape[1])//8, axis=1).mean().values
    beta532molAv = pd.DataFrame(beta532mol).groupby(np.arange(beta532mol.shape[1])//8, axis=1).mean().values
    const0 = np.mean(scc0Av[:,z_cc]/beta355molAv[:,z_cc], axis=1).reshape((time.size, 1))
    const1 = np.mean(scc1Av[:,z_cc]/beta532molAv[:,z_cc], axis=1).reshape((time.size, 1))
    scc0calib = scc0Av/const0
    scc1calib = scc1Av/const1

    N = np.int_(np.linspace(0, len(time)-1, 16))
    f0, axs = plt.subplots(ncols=4, nrows=4, figsize=(16,10))
    for n,ax in enumerate(axs.flatten()):
        # ax.plot(scc0calib[N[n],:], alt[::8])
        # ax.plot(beta355molAv[N[n],:], alt[::8], '--', color='red')
        ax.plot(scc0calib[N[n],:]/beta355molAv[N[n],:], alt[::8], label='355', color='blue')
        ax.plot(scc1calib[N[n],:]/beta532molAv[N[n],:], alt[::8], label='532', color='green')
        ax.vlines(1, ymin=alt[0], ymax=alt[-1], linestyles='--', color="black", zorder=10)
        ax.legend()
        ax.set(title=str(time[N[n]]))

    plt.suptitle(f'Scattering Ratio signal \n {fl}')
    plt.tight_layout()
    plt.savefig('/homedata/nmpnguyen/IPRAL/SCC_produits/sr'+fl.name.split('.')[0]+'.png')

    f1, axs = plt.subplots(ncols=4, nrows=4, figsize=(16,10))
    for n,ax in enumerate(axs.flatten()):
        ax.plot(rcs0[N[n],:], alt, label = '355:rcs')
        # ax.plot(scc0[N[n],:], alt, label = '355:rcs-atm')
        ax.legend()
        ax.set(title=str(time[N[n]]))

    plt.suptitle(f'Range_corrected_signal - atmospheric_background \n {list_scc[0]}')
    plt.tight_layout()
    plt.savefig('/homedata/nmpnguyen/IPRAL/SCC_produits/rcs-atm'+fl.name.split('.')[0]+'.png')

    print('Save as NetCDF')
    ds = xr.Dataset(
        {'atb_532': (('time', 'alt'), scc1calib),
        'atb_355': (('time', 'alt'), scc0calib),
        'atbmol_532': (('time','alt'), beta532molAv),
        'atbmol_355': (('time', 'alt'), beta355molAv),
        'rcs_532': (('time', 'alt'), scc1Av),
        'rcs_355': (('time', 'alt'), scc0Av)
        },
        coords={
            'time': time, 'alt': alt[::8]},
    )
    ds.to_netcdf(Path('/homedata/nmpnguyen/IPRAL/SCC_produits', fl.name.split('.')[0]+'_atb.nc'))

