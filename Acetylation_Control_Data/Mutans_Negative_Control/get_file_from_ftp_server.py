#!/usr/bin/env python3
# encoding: utf-8

import ursgal
import os


def main():
    """
    Downloads an example file from PRIDE

    usage:
        ./get_file_from_ftp_server_example.py


    This is an simple example which shows how to load data from a ftp server.


    """
    params = {
        "ftp_url": "ftp://ftp.pride.ebi.ac.uk",
        "ftp_folder": "pride/data/archive/2022/11/PXD038045",
        "ftp_include_ext": [
            "FE279LPMa_PG_Slot2-23_1_4159.d.raw.zip",
        ],
        "ftp_output_folder": "dataFolder",
    }
    if os.path.exists(params["ftp_output_folder"]) is False:
        print("Ceated folder: {0}".format(params["ftp_output_folder"]))
        os.mkdir(params["ftp_output_folder"])
    R = ursgal.UController(params=params)
    R.fetch_file(engine="get_ftp_files_1_0_0")
    return


if __name__ == "__main__":
    main()
