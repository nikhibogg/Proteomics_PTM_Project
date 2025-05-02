#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os
import pandas as pd

def main(folder=None, enzyme=None, target_decoy_database=None):
    """
    An example test script to search all mzML files which are present in the
    specified folder. The search is currently performed with 3 search engines.

    The digestion enzyme has to be specified as well as the target-decoy
    database.

    usage:

        ./protein_database_serach.py <mzML_folder> <enzyme> <target_decoy_database>

    """
    # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    for mzml in glob.glob(os.path.join("{0}".format(folder), "*.mzML")):
        mzML_files.append(mzml)

    # We specify all search engines and validation engines that we want to use in a list
    # (version numbers might differ on windows or mac):
    search_engines = [
        "msfragger_3_0",
        # "xtandem_vengeance",
        # "msgfplus_v2019_07_03",
    ]

    validation_engine = "percolator_3_4_0"

    # Modifications that should be included in the search
    all_mods = [
        "S,opt,any, Phospho",
        "T,opt,any, Phospho",
        "Y,opt,any, Phospho",
    ]

    # Initializing the Ursgal UController class with
    # our specified modifications and mass spectrometer
    params = {
        "database": target_decoy_database,
        "enzyme": enzyme,
        "modifications": all_mods,
        "precursor_mass_tolerance_minus": 5,
        "precursor_mass_tolerance_plus": 5,
        "precursor_mass_tolerance_unit": "ppm",
        "frag_mass_tolerance": 20,
        "frag_mass_tolerance_unit": "ppm",
        "csv_filter_rules": [
            ["Is decoy", "equals", "false"],
            ["PEP", "lte", 0.01],
        ],
        "cpus": 1,
        "-xmx": "4g",
    }

    uc = ursgal.UController(params=params)

    # complete workflow:
    # every spectrum file is searched with every search engine,
    # results are validated (for each engine seperately),
    # validated results are merged and filtered for targets and PEP <= 0.01.
    # In the end, all filtered results from all spectrum files are merged
    result_files = []
    for spec_file in mzML_files:
        validated_results = []
        for search_engine in search_engines:
            unified_search_results = uc.search(
                input_file=spec_file,
                engine=search_engine,
            )
            validated_csv = uc.validate(
                input_file=unified_search_results,
                engine=validation_engine,
            )
            validated_results.append(validated_csv)

        validated_results_from_all_engines = uc.execute_misc_engine(
            input_file=validated_results,
            engine="merge_csvs_1_0_0",
        )
        filtered_validated_results = uc.execute_misc_engine(
            input_file=validated_results_from_all_engines,
            engine="filter_csv_1_0_0",
        )
        result_files.append(filtered_validated_results)

    results_all_files = uc.execute_misc_engine(
        input_file=result_files,
        engine="merge_csvs_1_0_0",
    )
    print(result_files)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        sys.exit(1)
    main(
        folder=sys.argv[1],
        enzyme=sys.argv[2],
        target_decoy_database=sys.argv[3],
    )

# How effective the tools are based on the type of tools they are.
# Which predictions are also verified, which families of proteins are more liekly to be identified.
# functional categories of proteins.
# Look at where sides that are correctly predicted are.
# Ratios of T S and Y identified.
# COmpare tool output, which tools based on the type 