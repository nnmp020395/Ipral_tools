#!/usr/bin/env python3
"""Script to convert IPRAL netCDF 1A with background corrected signal."""

import os
from pathlib import Path
import shutil
import sys

import click
import netCDF4 as nc
import numpy as np

__author__ = "marc-antoine drouin"
__email__ = "marc-antoine.drouin@lmd.ipsl.fr"
__version__ = "0.1.0"


IPRAL_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/ipral/")
IPRAL_MASK = "ipral_1a_Lz1R15mF30sPbck_v01_*_000000_1440.nc"


def substract_bckgrd(rcs, bckgrd, r_square):
    """
    Remove background from Pr2.

    Parameters
    ----------
    rcs : array(m, n)
        RCS profiles data.
    bckgrd : array(m)
        Background signal array.
    r_square : array(n)
        Range power 2.

    Returns
    -------
    np.array(m, n)
        RCS profiles with background substracted.
    """
    data = ((rcs / r_square).T - bckgrd).T * r_square

    return data


# @click.command()
# @click.argument("date", type=click.DateTime(formats=["%Y-%m-%d", "%Y%m%d"]))
# @click.argument(
#     "out_dir", type=click.Path(dir_okay=True, writable=True, resolve_path=True), nargs=1
# )
def convert_rcs(data_file, date, out_dir):
    """
    Convert in IPRAL 1a file RCS variables to background corrected signal.

    Parameters
    ----------
    date : datetime.datetime
        Date to process.
    out_dir : str or pathlib.Path
        Directory where the file will be store.

    """
    # out_dir = Path(out_dir)
    # data_dir = IPRAL_PATH / f"{date:%Y}" / f"{date:%m}" / f"{date:%d}"
    # data_file = sorted(data_dir.glob(IPRAL_MASK))[0]

    # try to copy file
    print(f"copying {data_file.name} into {out_dir}")
    try:
        new_file = shutil.copy(data_file, out_dir)
    except FileNotFoundError:
        print(f"ERROR: no IPRAL file available for {date:%Y-%m-%d}")
        return 1

    # change right on the file
    print(f"change right of file {new_file}")
    os.chmod(new_file, 0o666)

    # change RCS values of RCS variables
    nc_id = nc.Dataset(new_file, "a")

    r_square = np.square(nc_id.variables["range"][:])

    for nc_var in nc_id.variables:
        # we only search for variables starting with RCS
        if not nc_var.startswith("rcs_"):
            continue

        print(f"substracting background to variable {nc_var}")
        rcs = nc_id.variables[nc_var][:]
        index = nc_var.split("_")[-1]
        bckgrd = nc_id.variables["bckgrd_rcs_" + index][:]
        # print(f'print bckgrd-----------')
        # print(bckgrd)
        nc_id.variables[nc_var][:] = substract_bckgrd(rcs, bckgrd, r_square)
        nc_id.variables[
            nc_var
        ].long_name = "range corrected signal background substracted"

    print("end")

    return 0


if __name__ == "__main__":
    sys.exit(convert_rcs())
