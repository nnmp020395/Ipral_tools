import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from pathlib import Path
import sys
import xarray as xr

import numpy as np
import pandas as pd

simul_file = Path("/homedata/nmpnguyen/IPRAL/1a/ipral_1a_Lz1R15mF30sPbck_v01_20200914_000000_1440_simul.pkl")
ipral_file = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/2020/09/14/ipral_1a_Lz1R15mF30sPbck_v01_20200914_000000_1440.nc")
# "/homedata/nmpnguyen/IPRAL/1a/ipral_1a_bck_corrected_v01_20200121_000000_1440.nc"
# sftp://nmpnguyen@camelot.ipsl.polytechnique.fr/mnt/sirtafs-2.5/lov/
print(simul_file, ipral_file)
def get_calib_plot(rcs17atb, beta532molAv, rcs13atb, beta355molAv, time, rAv, ipral_file, wave):
    N = np.int_(np.linspace(0, len(time)-1, 16))
    rplot = rAv[np.where(rAv<=20000)]
    fg, axs = plt.subplots(figsize=(16,12), nrows=4, ncols=4, sharey=True)
    if wave == 355:
        for m, ax in enumerate(axs.flatten()):
            sr = rcs13atb.iloc[N[m],:rplot.shape[0]]/beta355molAv.iloc[N[m],:rplot.shape[0]]
            ax.plot(sr, rplot, label="SR 355nm")
            ax.vlines(1, ymin=rplot[0], ymax=rplot[-1], linestyles='--', color="red", zorder=10)
            # ax.plot(rcs13atb.iloc[N[m],:rplot.shape[0]], rplot, label="ATB 355nm")
            # ax.plot(beta355molAv.iloc[N[m],:rplot.shape[0]], rplot, "--", label="mol 355nm", color="red")
            ax.axhline(3900, linestyle='--', color='black')   
            ax.set(ylabel="Alt,m", title=str(time[N[m]]))
            ax.legend()
            # ax.set_ylim(0, 20000)
        # plt.gca().set_xlim(left=0)
        plt.suptitle(str(ipral_file.name)+'\n Attenuated Backscatter Profile')
        plt.tight_layout()        
        out_png = Path('/homedata/nmpnguyen/IPRAL/1a_Fig', ipral_file.stem + '_355sr.png')
        plt.savefig(out_png)
    else:
        for m, ax in enumerate(axs.flatten()):
            sr = rcs17atb.iloc[N[m],:rplot.shape[0]].values/beta532molAv.iloc[N[m],:rplot.shape[0]].values
            ax.plot(sr, rplot, label="SR 532nm")
            ax.vlines(1, ymin=rplot[0], ymax=rplot[-1], linestyles='--', color="red", zorder=10)
            # ax.plot(rcs17atb.iloc[N[m],:rplot.shape[0]], rplot, label="ATB 532nm")
            # ax.plot(beta532molAv.iloc[N[m],:rplot.shape[0]], rplot, "--", label="mol 532nm", color="red")
            ax.axhline(3900, linestyle='--', color='black')            
            ax.set(ylabel="Alt,m", title=str(time[N[m]]))
            ax.legend()
            # ax.set_ylim(0, 20000)
        # plt.gca().set_xlim(left=0)
        plt.suptitle(str(ipral_file.name)+'\n Scattering Ratio profiles')
        plt.tight_layout()
        out_png = Path('/homedata/nmpnguyen/IPRAL/1a_Fig', ipral_file.stem + '_532sr.png')
        plt.savefig(out_png)
        # plot quicklook atb    
        fig, ax = plt.subplots()
        p = ax.pcolormesh(time, rAv, rcs17atb.T, norm=colors.LogNorm(vmin=10e-9, vmax=10e-6), cmap = 'RdBu_r', shading='nearest')
        plt.colorbar(p, ax=ax)
        ax.set_ylim(0,20000)
        plt.suptitle(str(ipral_file.name)+'\n Attenuated Backscatter Quicklook')
        plt.tight_layout()
        out_png = Path('/homedata/nmpnguyen/IPRAL/1a_Fig', ipral_file.stem + '_532atbQL.png')
        plt.savefig(out_png)

"""
Fontionc init et add sont pour détecter le pic du signal, donc détecter les nuages, aérosols
"""

def init(
    x,
    lag,
    threshold,
    influence,
    ):
    '''
    Smoothed z-score algorithm
    Implementation of algorithm from https://stackoverflow.com/a/22640362/6029703
    '''
    labels = np.zeros(lag)
    filtered_y = np.array(x[0:lag])
    avg_filter = np.zeros(lag)
    std_filter = np.zeros(lag)
    var_filter = np.zeros(lag)
    avg_filter[lag - 1] = np.mean(x[0:lag])
    std_filter[lag - 1] = np.std(x[0:lag])
    var_filter[lag - 1] = np.var(x[0:lag])
    return dict(avg=avg_filter[lag - 1], var=var_filter[lag - 1],
                std=std_filter[lag - 1], filtered_y=filtered_y,
                labels=labels)


def add(result, 
    single_value,
    lag,
    threshold,
    influence,
    ):
    previous_avg = result['avg']
    previous_var = result['var']
    previous_std = result['std']
    filtered_y = result['filtered_y']
    labels = result['labels']
    if abs(single_value - previous_avg) > threshold * previous_std:
        if single_value > previous_avg:
            labels = np.append(labels, 1)
        else:
            labels = np.append(labels, -1)
        # calculate the new filtered element using the influence factor
        filtered_y = np.append(filtered_y, influence * single_value
                               + (1 - influence) * filtered_y[-1])
    else:
        labels = np.append(labels, 0)
        filtered_y = np.append(filtered_y, single_value)
    # update avg as sum of the previuos avg + the lag * (the new calculated item - calculated item at position (i - lag))
    current_avg_filter = previous_avg + 1. / lag * (filtered_y[-1]
            - filtered_y[len(filtered_y) - lag - 1])
    # update variance as the previuos element variance + 1 / lag * new recalculated item - the previous avg -
    current_var_filter = previous_var + 1. / lag * ((filtered_y[-1]
            - previous_avg) ** 2 - (filtered_y[len(filtered_y) - 1
            - lag] - previous_avg) ** 2 - (filtered_y[-1]
            - filtered_y[len(filtered_y) - 1 - lag]) ** 2 / lag)  # the recalculated element at pos (lag) - avg of the previuos - new recalculated element - recalculated element at lag pos ....
    # calculate standard deviation for current element as sqrt (current variance)
    current_std_filter = np.sqrt(current_var_filter)
    return dict(avg=current_avg_filter, var=current_var_filter,
                std=current_std_filter, filtered_y=filtered_y[1:],
                labels=labels)


def get_calibration(ipral_file, simul_df):
    # from timeit import default_timer as timer
    # TIME = timer()
    # SIMULATED DATA------------------------------- input simulated variables from era5
    # new_df = pd.read_pickle(simul_file)
    beta532mol = simul_df['beta532mol']
    beta532mol = beta532mol.unstack(level=1)
    beta355mol = simul_df['beta355mol']
    beta355mol = beta355mol.unstack(level=1)
    print('SIMULATED DATA-------------------------------end')
    # IPRAL DATA-------------------------------
    d = xr.open_dataset(ipral_file)
    time = d.time.values
    r = d.range.values
    r_square = np.square(r)
    ipr17 = d["rcs_17"].values
    bck17 = d['bckgrd_rcs_17'].values
    tmp = ipr17/r_square.reshape((1, r_square.size)) - bck17.reshape((bck17.size, 1))
    rcs17 = tmp*r_square.reshape((1, r_square.size))
    rcs17 = pd.DataFrame(data=rcs17, index=time, columns=np.float64(r))
    ipr13 = d["rcs_13"].values
    bck13 = d['bckgrd_rcs_13'].values
    tmp = ipr13/r_square.reshape((1, r_square.size)) - bck13.reshape((bck13.size, 1))
    rcs13 = tmp*r_square.reshape((1, r_square.size))
    rcs13 = pd.DataFrame(data=rcs13, index=time, columns=np.float64(r))
    # print(f'Timing : {timer()-TIME}')
    print('IPRAL RANGE CORRECTED DATA-------------------------------end')
    # AVERAGE-------------------------------
    rcs17av = rcs17.groupby(np.arange(len(rcs17.columns))//8, axis=1).mean()
    rcs13av = rcs13.groupby(np.arange(len(rcs13.columns))//8, axis=1).mean()
    beta532molAv = beta532mol.groupby(np.arange(len(beta532mol.columns))//8, axis=1).mean()
    beta355molAv = beta355mol.groupby(np.arange(len(beta355mol.columns))//8, axis=1).mean()
    rcs17av.columns = r[::8]
    rcs13av.columns = r[::8]
    beta532molAv.columns = r[::8]
    beta355molAv.columns = r[::8]
    rAv = r[::8]
    print('AVERAGE-------------------------------end')
    # PEAK DETECTION
    r_id = np.where(rAv <= 4000)[0]
    time_mask = []
    id_mask = []
    for t in range(time.size):
        profil1 = rcs17av.iloc[t,r_id].values
        lag = 5
        threshold = 5
        influence = 0.5
        m=0
        result = init(profil1, lag=lag, threshold=threshold, influence=influence)
        for line in profil1:
            result = add(result, line, lag, threshold, influence)
        if (len(np.where(result['labels']==1)[0]) != 0):
            time_mask.append(time[t])
            id_mask.append(t)
    #---------        
    print(f'time index of peak detected: {id_mask}')
    rcs17av_new = np.array(rcs17av, copy=True)
    rcs17av_new[id_mask,:] = np.nan
    rcs17av_new = pd.DataFrame(rcs17av_new, index = time, columns=rAv)
    print('PEAK DETECTION-------------------------end')
    # NORMALIZATION-------------------------
    z_cc = np.where((rAv>3700)&(rAv<4000))[0]
    constk = rcs17av_new.iloc[:,z_cc].div(beta532molAv.iloc[:,z_cc]).mean(axis=1) 
    rcs17atb = rcs17av_new.div(constk, axis=0)
    constk = rcs13av.iloc[:,z_cc].div(beta355molAv.iloc[:,z_cc]).mean(axis=1)
    rcs13atb = rcs13av.div(constk, axis=0)
    sr17 = rcs17atb.div(beta532molAv)
    sr13 = rcs13atb.div(beta355molAv)
    print('NORMALIZATION-------------------------end')
    get_calib_plot(rcs17atb, beta532molAv, rcs13atb, beta355molAv, time, rAv, ipral_file, 355)
    get_calib_plot(rcs17atb, beta532molAv, rcs13atb, beta355molAv, time, rAv, ipral_file, 532)
    ds = xr.Dataset(
        {'atb_17': (('time', 'range'), rcs17atb.values),
        'atb_13': (('time', 'range'), rcs13atb.values),
        'atbmol_17': (('time','range'), beta532molAv.values),
        'atbmol_13': (('time', 'range'), beta355molAv.values),
        'rcs_17': (('time', 'range'), rcs17av),
        'rcs_13': (('time', 'range'), rcs13av)},
        coords={
            'time': time, 'range': rAv},
    )
    ds.to_netcdf(Path('/homedata/nmpnguyen/IPRAL/1a', ipral_file.name.split('.')[0]+'_atb.nc'))
    return rAv, rcs17atb, rcs13atb, sr17, sr13, beta355molAv, beta532molAv, time
    # return r, rcs17atb, rcs13atb, sr17, sr13, beta355mol, beta532mol, time


get_calibration(ipral_file, pd.read_pickle(simul_file))


# ax[0].vlines(1, ymin=0, ymax=20000, linestyles="--", color="red", zorder=10)
# ax[0].plot(rcs17atb.iloc[m,:], r, label="ATB 532nm")
# ax[0].plot(beta532mol.iloc[m,:], r, "--", label="mol 532nm", color="red")
# # ax[0].plot(rcs13atb.iloc[m,:], r, label="ATB 355nm", color="green")
# # ax[0].plot(beta355mol.iloc[m,:], r, "--", label="mol 355nm", color="black")
# ax[0].set(ylabel="Alt,m", xlabel="ATB")
# ax[0].legend()
# ax[0].set_ylim(0, 20000)
# ax[1].plot(sr17.iloc[m,:], r, label="SR 532nm")
# # ax[1].plot(sr13.iloc[m,:], r, label="SR 355nm", color="green")
# ax[1].vlines(1, ymin=0, ymax=20000, linestyles="--", color="red", zorder=10)
# ax[1].set(xlabel="SR")
# ax[1].legend()
# ax[1].set_ylim(0, 20000)
# ax[1].set_xlim(0,3)
# plt.suptitle(str(ipral_file.name) + "\n" + str(time[m]))
# plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/rcs13profiles.png")
