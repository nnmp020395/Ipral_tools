#!/usr/bin/env python3
"""Script to ananlyze IPRAL L1a files."""

import datetime as dt
from pathlib import Path
import sys

import click
import xarray as xr

__author__ = "marc-antoine drouin"
__email__ = "marc-antoine.drouin@lmd.ipsl.fr"
__version__ = "0.1.0"


IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"


@click.command()
@click.argument("date_start", type=click.DateTime(formats=["%Y-%m-%d", "%Y%m%d"]))
@click.argument("date_end", type=click.DateTime(formats=["%Y-%m-%d", "%Y%m%d"]))
def analyze_ipral_file(date_start, date_end):
    """
    Analyze IPRAL L1a file.

    The script returns filename, first and last timestep and number of profiles.

    Parameters
    ----------
    date_start : datetime.datetime
        First date to search.
    date_end : datetime.datetime
        Last date to search.

    """
    # list of dates to search files in
    list_dates = [
        date_start + dt.timedelta(days=d)
        for d in range((date_end - date_start).days + 1)
    ]

    for date in list_dates:

        files_path = IPRAL_PATH / f"{date:%Y}" / f"{date:%m}" / f"{date:%d}"
        list_files = sorted(files_path.glob(IPRAL_MASK))

        if not list_files:
            continue

        time = (
            xr.open_dataset(list_files[0])["time"].to_dataframe().index.to_pydatetime()
        )

        print(
            f"{list_files[0].name} {time[0]:%Y-%m-%dT%H:%M:%SZ} {time[-1]:%Y-%m-%dT%H:%M:%SZ} {time.size:4d}"  # NOQA
        )

    return 0


if __name__ == "__main__":
    sys.exit(analyze_ipral_file())
