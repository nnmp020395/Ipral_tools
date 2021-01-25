#!/usr/bin/env python3
"""Script to ananlyze IPRAL L1a files"""

import datetime as dt
from pathlib import Path
import sys

import click
import netCDF4 as nc

__author__ = "marc-antoine drouin"
__email__ = "marc-antoine.drouin@lmd.ipsl.fr"
__version__ = "0.1.0"


IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"


class DateParamType(click.ParamType):
    """Implement a date type param for click package."""

    name = "date"

    def __init__(self, format):
        """Add format of date."""
        self.format = format

    def convert(self, value, param, ctx):
        """Convert date string into datetime object of possible."""
        try:
            date_dt = dt.datetime.strptime(value, self.format)
        except ValueError:
            self.fail(
                "%s is not a valid date with format %s" % (value, self.format),
                param,
                ctx,
            )

        return date_dt


@click.command()
@click.argument("date_start", type=DateParamType("%Y-%m-%d"), nargs=1)
@click.argument("date_end", type=DateParamType("%Y-%m-%d"), nargs=1)
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

        nc_id = nc.Dataset(list_files[0])
        time = nc.num2date(
            nc_id.variables["time"][:],
            units=nc_id.variables["time"].units,
            only_use_cftime_datetimes=False,
        )
        nc_id.close()

        print(
            f"{list_files[0].name} {time[0]:%Y-%m-%dT%H:%M:%SZ} {time[-1]:%Y-%m-%dT%H:%M:%SZ} {time.size:4d}"  # NOQA
        )

    return 0


if __name__ == "__main__":
    sys.exit(analyze_ipral_file())
