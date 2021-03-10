import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import re
from pathlib import Path
import sys
import xarray as xr
import numpy as np
import pandas as pd

sys.path.append('/homedata/nmpnguyen/ipral-tools/')
from ipral_1a_bck_corrected import convert_rcs
from ipral_chm15k_cloud_filter import ipral_remove_cloud_profiles
import datetime as dt
# from imp import reload as rl
# rl(ipral_variables_simulation)
from ipral_variables_simulation import simulate 
from calibrer_test import get_calibration, get_calib_plot
#(1): Find all ipral file in folder 

"""
background correction ==> convert_rcs(ipral_file, date, outdir)
cloud_filter ==> ipral_remove_cloud_profiles(date, alt_max, ipral_file, outdir)
"""
IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"

input_year = '2020'
data_dir = IPRAL_PATH / input_year
list_files = sorted(data_dir.rglob(IPRAL_MASK))

for lf in list_files[:5]:# '/bdd/SIRTA/pub/basesirta/1a/ipral/2019/05/14/ipral_1a_Lz1R15mF30sPbck_v01_20190514_000000_1440.nc'
# for lf, ax in zip(list_files[16:], axs.flatten()):
    print(lf)
    outdir_parent = "/homedata/nmpnguyen/IPRAL/1a"
    date = dt.datetime.strptime(lf.name.split('_')[4], '%Y%m%d')  
    out_filter = outdir_parent / Path(re.sub(lf.name.split('_')[2], 'cloud_filter', lf.name))
    ipral_remove_cloud_profiles(date, 3700, lf, out_filter)
    # out_bck = Path(re.sub('cloud_filter','bck_corrected',str(out_filter))) #outdir_parent / Path('_'.join(index)) 
    # convert_rcs(lf, date, out_bck)
    simul_df = simulate(out_filter, 'era')
    rAv, rcs17atb, rcs13atb, sr17, sr13, beta355molAv, beta532molAv, time=get_calibration(out_filter, simul_df)
    
    # get_calib_plot(rcs17atb, beta532molAv, rcs13atb, beta355molAv, time, rAv, out_filter, 532)
    # get_calib_plot(rcs17atb, beta532molAv, rcs13atb, beta355molAv, time, rAv, out_filter, 355)
#     mol355_mean_day = beta355molAv.mean(axis=0)
#     atb355_mean_day = rcs13atb.mean(axis=0)
#     zr355_mean_day = atb355_mean_day.div(mol355_mean_day)
#     mol532_mean_day = beta532molAv.mean(axis=0)
#     atb532_mean_day = rcs17atb.mean(axis=0)
#     zr532_mean_day = atb532_mean_day.div(mol532_mean_day)  
#     ax.plot(zr532_mean_day, rAv, label="532", color="blue")
#     # ax.plot(zr355_mean_day, rAv, label="355", color="green")
#     ax.vlines(1, ymin=rAv[0], ymax=rAv[-1], linestyles='--', color="red", zorder=10)
#     ax.legend(loc="upper center")
#     ax.set(title=str(date))
#     ax.set_ylim(0, 22000)
#     ax.set_xlim(0, 5)

# plt.suptitle(f'Scattering Ratio signal')
# plt.tight_layout()
# plt.savefig('/homedata/nmpnguyen/IPRAL/1a_Fig/Ipral_SR_average.png')