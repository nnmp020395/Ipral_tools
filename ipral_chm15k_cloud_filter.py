#!/usr/bin/env python3
"""
Script to remove IPRAL profiles with clouds.

The script use data from SIRTA CHM15k and a maximum altitude
of clouds defined by the user.
"""
from pathlib import Path
import sys

import click
import xarray as xr


__author__ = "marc-antoine drouin"
__email__ = "marc-antoine.drouin@lmd.ipsl.fr"
__version__ = "1.0.0"


CHM15K_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/chm15k")
CHM15K_MASK = "chm15k_1a_z1Ppr2R15mF15s_v01_*_1440.nc"
CHM15K_TIME_RES = "15s"

IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_1440.nc"
IPRAL_TIME_RES = "30s"

DATE_FMT = "days since 1970-01-01T00:00:00"


def find_file(sirta_data_path, data_mask, date):
    """
    Search files in SIRTA database from root path and mask.

    Parameters
    ----------
    sirta_data_path : pathlib.Path
        The root path of the data in SIRTA database.
    data_mask : str
        Mask of the fileanme to search for.
    date : datetime.datetime
        The date of the file to search.

    Returns
    -------
    list of pathlib.Path
        The sorted list of files found.

    """
    data_path = sirta_data_path / f"{date:%Y}" / f"{date:%m}" / f"{date:%d}"

    data_files = sorted(data_path.glob(data_mask))

    return data_files


@click.command()
@click.argument("date", type=click.DateTime(formats=["%Y-%m-%d", "%Y%m%d"]))
@click.argument("alt-max", type=click.FLOAT)
@click.argument(
    "ipral_file",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "output",
    type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True),
)
def ipral_remove_cloud_profiles(date, alt_max, ipral_file, output):
    """
    Remove IPRAL profiles containing cloud below a defined altitude.

    Parameters
    ----------
    date : datetime.datetime
        The date of the file to process.
    alt_max : float
        The maximum altitude of clouds in meters.
    ipral_file : str or pathlib.Path
        The path of the IPRAL file to process.
    output : str or pathlib.Path
        The path to th eoutput file.

    """
    print(f"Processing {date:%Y-%m%d}")
    print(f"Removing IPRAL profiles with clouds below {alt_max:7.1f}")
    # read CHM15k file
    # ---------------------------------------------------------------------------------
    chm15k_file = find_file(CHM15K_PATH, CHM15K_MASK, date)
    if not chm15k_file:
        print("No CHM15k file found.")
        print("Quitting.")
        sys.exit(1)

    chm15k_file = chm15k_file[0]
    print(f"CHM15k file found: {chm15k_file}")

    cbh = xr.open_dataset(chm15k_file)["cloud_base_height"][:, 0].to_dataframe()[
        "cloud_base_height"
    ]
    # round time to 15s to ease use
    cbh.index = cbh.index.round(freq=CHM15K_TIME_RES)
    # under sample chm15k data to 30s to have the time resolution as ipral
    cbh = cbh.resample(IPRAL_TIME_RES).first()

    # read IPRAL data
    # ---------------------------------------------------------------------------------
    ipral_data = xr.open_dataset(ipral_file)
    raw_profs = ipral_data.time.size
    print(f"{raw_profs} in IPRAL data")

    # get cloud mask
    # ---------------------------------------------------------------------------------
    # round time to 30s to ease use
    ipral_time = ipral_data.time.to_dataframe().index.round(freq=IPRAL_TIME_RES)
    # only keep timesteps of CBH available in ipral data
    cbh = cbh.loc[ipral_time]
    # create to only keep data without cloud under the chosen altitude
    cbh_mask = cbh > alt_max
    profs_to_keep = cbh_mask.values.astype("i2").sum()
    print(f"{raw_profs - profs_to_keep} profiles will be remove")

    # apply mask
    # ---------------------------------------------------------------------------------
    ipral_data = ipral_data.isel(time=cbh_mask)

    # save file
    # ---------------------------------------------------------------------------------
    print(f"saving in {output}")
    ipral_data.to_netcdf(
        output,
        format="NETCDF4",
        encoding={"time": {"units": DATE_FMT, "calendar": "standard"}},
    )

    return 0


if __name__ == "__main__":
    sys.exit(ipral_remove_cloud_profiles())
