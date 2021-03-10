import matplotlib
matplotlib.use('Agg')
import sys 
# sys.path.append("/homedata/nmpnguyen/ipral-tools/")
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
import click 

def trim_axs(axs, N):
    """
    Reduce *axs* to *N* Axes. All further Axes are removed from the figure.
    """
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]
# import ipral_1a_bck_corrected
#input file
@click.command()
# @click.argument("file_ipral", type=click.Path(dir_okay=True, writable=True, resolve_path=True))
@click.argument("date", type=click.DateTime(formats=["%Y-%m-%d", "%Y%m%d"]))

# fil((e_ipral = Path("/homedata/nmpnguyen/IPRAL/1a/2019-05-14_ipral_1a_bck_corrected.nc")
# print(date)
def plot_raw(date):
    print(f'Plot Date file: {date}')
    # load data------------------
    # IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
    # IPRAL_DIR = IPRAL_PATH / f"{date:%Y}" / f"{date:%m}" / f"{date:%d}"
    # IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"
    IPRAL_DIR = Path('/homedata/nmpnguyen/IPRAL/1a/')
    IPRAL_MASK = "ipral_1a_cloud_filter_v01_"+date.strftime('%Y%m%d')+"_000000_1440.nc"
    file_ipral = sorted(IPRAL_DIR.glob(IPRAL_MASK))[0]
    print(f'Datafile path: {file_ipral}')
    # file_ipral = Path("/homedata/nmpnguyen/IPRAL/1a")/Path(date+"_ipral_1a_bck_corrected.nc")
    data = xr.open_dataset(file_ipral)
    # plot rcs_13----------------
    fig, ax = plt.subplots(figsize=[10, 6])
    np.log(data.rcs_13).plot(ax=ax, x='time', vmin=15, vmax=20)
    # plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/"+date+"_ipral_1a_bck_corrected_rcs13.png")
    ax.set_ylim(0, 20000)
    plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/"+str(date)+"_ipral_1a__rcs13.png")
    plt.close(fig)
    # plot rcs_17----------------
    fig, ax = plt.subplots(figsize=[10, 6])
    np.log(data.rcs_17).plot(ax=ax, x='time', vmin=15, vmax=20)
    # plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/"+date+"_ipral_1a_bck_corrected_rcs17.png")
    ax.set_ylim(0, 20000)
    plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/"+str(date)+"_ipral_1a__rcs17.png")
    plt.close(fig)

    
# # Profiles plot rcs_13----------------
def profiles(date, channel, file_path):
    if file_path == None : 
        sys.path.append('/homedata/nmpnguyen/ipral-tools/')
        from ipral_1a_bck_corrected import convert_rcs
        IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
        IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"
        date = dt.datetime.strptime(date, '%Y-%m-%d')
        try:
            ipral_file = sorted(Path(IPRAL_PATH, f'{date:%Y}', f"{date:%m}", f"{date:%d}").glob(IPRAL_MASK))[0]
        except IndexError:
            print(f"ERROR: no IPRAL file available for {date:%Y-%m-%d}")
            return 1
    else: 
        ipral_file = Path(file_path)
    #----   
    print(f'Open IPRAL file')
    d = xr.open_dataset(ipral_file)
    time = d.time.values
    r = d.range.values
    r_square = np.square(r)
    channel_rcs = 'rcs_'+str(channel)
    channel_bckgrd = 'bckgrd_rcs_'+str(channel)
    ipr = d[channel_rcs].values
    bck = d[channel_bckgrd].values
    tmp = ipr/r_square.reshape((1, r_square.size)) - bck.reshape((bck.size, 1))
    rcs = tmp*r_square.reshape((1, r_square.size))
    rcs = pd.DataFrame(data=rcs, index=time, columns=np.float64(r))
    #----
    fig, axs = plt.subplots(figsize=(16,12), ncols=5, nrows=5, sharey=True, sharex=True)
    N = np.int_(np.linspace(0, len(time)-1, 25))
    for n, ax in enumerate(axs.flatten()):
        ax.plot(ipr[N[n],:], r, label=channel_rcs)
        ax.plot(rcs.loc[time[N[n]]], r, label=channel_bckgrd)
        ax.set(title=str(time[N[n]]))
        ax.legend()
        ax.set_ylim(0, 20000)
    plt.suptitle(str(ipral_file.name) + "\n Range background correction")
    plt.tight_layout()
    plt.savefig("/homedata/nmpnguyen/IPRAL/1a_Fig/"+ipral_file.name.split('.')[0]+channel_rcs+'_profile.png')


if __name__ == "__main__":
    sys.exit(plot_raw())
