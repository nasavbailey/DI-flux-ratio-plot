# Generate Roman Coronagraph predicted performance curves.
# This script requires corgietc to be installed.
# See: https://github.com/roman-corgi/corgietc for installation instructions

import glob
import os
from astropy.io import ascii
import pandas
import corgietc  # noqa
import json
import copy
import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping
import astropy.units as u
from astropy.table import Table
import numpy as np

# identify all current curves
datapath = "./data/"
curr_curves = glob.glob(os.path.join(datapath, "Roman_pred*.txt"))

fnames = []
names = []
lams = []
tints = []
fpps = []
SNRs = []
for c in curr_curves:
    fname = os.path.split(c)[-1]
    fnames.append(fname)
    names.append(fname.split("Roman_pred_")[-1].split(".txt")[0])
    dat = ascii.read(c)
    lams.append(dat["lambda"][0])
    tints.append(dat["t_int_hr"][0])
    fpps.append(dat["fpp"][0])
    SNRs.append(dat["SNR"][0])

data = pandas.DataFrame(
    {
        "Filename": fnames,
        "Name": names,
        "lambda": lams,
        "t_int_hr": tints,
        "fpp": fpps,
        "SNR": SNRs,
    }
)
# backup of original curve data for historical purposes
# data.to_csv(os.path.join(datapath, "Roman_CDR_curves_list.csv"), index=False)

# set up objects
scriptfile = os.path.join(os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json")
with open(scriptfile, "r") as f:
    specs = json.loads(f.read())

TK = EXOSIMS.Prototypes.TimeKeeping.TimeKeeping(
    missionLife=5.25
)  # 63 months in years is 5.25, 21 months is 1.75
TK.allocate_time(21 * 30.4375 * u.d)
TL = EXOSIMS.Prototypes.TargetList.TargetList(**copy.deepcopy(specs))
OS = TL.OpticalSystem

# figure out equivalent observing modes
scenarios = []
modes = np.array([mode["Scenario"] for mode in OS.observingModes])
for _, row in data.iterrows():
    scenario = "CON_" if row.Name.endswith("_cons") else "OPT_"
    if row.Name.startswith("imaging_"):
        scenario += "IMG_NFB1_HLC"
    elif row.Name.startswith("spec_"):
        scenario += "SPEC_NFB3_SPC"
    elif row.Name.startswith("wideFOVimaging_"):
        scenario += "IMG_WFB4_SPC"

    # consistency check
    assert scenario in modes, f"Could not match scenario for {row.Name}"
    assert (
        OS.observingModes[np.where(modes == scenario)[0][0]]["lam"].to_value("nm")
        == row["lambda"]
    ), f"Incorrect wavelength for matched scenario for {row.Name}"

    scenarios.append(scenario)

data["Scenario"] = scenarios

# set reused values
sInds = 0
fZ = np.repeat(TL.ZodiacalLight.fZ0, 1)
mode = OS.observingModes[0]
JEZ = TL.JEZ0[mode["hex"]] / (4.1536**2)

# compute curves
for _, row in data.iterrows():
    mode = OS.observingModes[np.where(modes == row.Scenario)[0][0]]
    assert mode["Scenario"] == row.Scenario

    WAs = (
        np.linspace(mode["IWA"].value * 1.01, mode["OWA"].value * 0.99, 100)
        * mode["IWA"].unit
    )

    # treat 10k hours as saturation
    if row.t_int_hr == 10000:
        dMags = OS.calc_saturation_dMag(
            TL,
            [sInds] * len(WAs),
            np.repeat(fZ, len(WAs)),
            np.repeat(JEZ, len(WAs)),
            WAs,
            mode,
            TK=TK,
        )
    # otherwise, use the actual time
    else:
        dMags = OS.calc_dMag_per_intTime(
            np.ones(len(WAs)) * row.t_int_hr * u.hr,
            TL,
            [sInds] * len(WAs),
            np.repeat(fZ, len(WAs)),
            np.repeat(JEZ, len(WAs)),
            WAs,
            mode,
            TK=TK,
        )

    # generate output table
    df = pandas.DataFrame(
        {
            "l/D": (WAs / mode["syst"]["input_angle_unit_value"]).value,
            "contr": 10 ** (-0.4 * dMags),
            "lambda": np.ones(len(WAs)) * mode["lam"].to_value("nm"),
            "t_int_hr": np.ones(len(WAs)) * row.t_int_hr,
            "fpp": np.ones(len(WAs)) * mode["pp_Factor_CBE"],
            "SNR": np.ones(len(WAs)) * mode["SNR"],
        }
    )

    # write to disk
    Table.from_pandas(df).write(
        os.path.join(datapath, row.Filename), format="ascii.ecsv", overwrite=True
    )
