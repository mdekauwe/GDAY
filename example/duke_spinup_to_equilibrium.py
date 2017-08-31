#!/usr/bin/env python

""" Duke Simulation for NCEAS FACE experiment

Site History:
-------------
* Pre-1850 temperate broadleaved deciduous forest
* Harvest forest.
* Grassland from 1850-1982, mowed (+/- annually - DO NOT SIMULATE MOWING :P)
* burnt prior to planting
* At start of experiment aboveground biomass 5.5-11 kg C m-2
* Needleaf forest planted in 1983

Spin-up the model to a steady state. Recycle the met data in batches of a
50 years (over and over), until the SOM, plant and litter C pools cease to
change.

-> Spinup with forest params, fixed NDEP, fixed CO2 .
-> Vary NDEP/CO2 for about 200 odd yrs...using grassland params so that we get
   through the industrial to the 1980s period.
"""

import os
import shutil
import sys
import subprocess
import numpy as np

USER = os.getlogin()
sys.path.append('/Users/%s/src/c/gday/scripts' % (USER))
import adjust_gday_param_file as ad

__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.12.2014)"
__email__   = "mdekauwe@gmail.com"


def main(experiment_id, site, SPIN_UP=True, POST_INDUST=True, SPIN_UP_SIMS=True):

    GDAY_SPIN = "gday -s -p "
    GDAY = "gday -p "

    base_dir = os.path.dirname(os.getcwd())

    # dir names
    base_dir = "/Users/%s/src/c/gday/example/params" % (USER)
    param_dir = "params"
    met_dir = os.path.join(base_dir, "met_data")
    run_dir = os.path.join(base_dir, "outputs")


    if SPIN_UP == True:

        # copy base files to make two new experiment files
        shutil.copy(os.path.join(base_param_dir, base_param_name + ".cfg"),
                    os.path.join(param_dir, "%s_%s_model_spinup.cfg" % \
                    (experiment_id, site)))

        # Run model to equilibrium assuming forest, growing C pools from effectively
        # zero
        itag = "%s_%s_model_spinup" % (experiment_id, site)
        otag = "%s_%s_model_spunup" % (experiment_id, site)
        mtag = "%s_met_data_equilibrium_50_yrs.csv" % (site)
        out_fn = itag + "_equilib.out"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {

                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),

                         # state - default C:N 25.
                         "age": "0.0",
                         "canht": "17.0", # Canopy height increased from 16m in 2001 to 18m in 2004 at Duke
                         "activesoil": "0.001",
                         "activesoiln": "0.00004",
                         "age": "0.0",
                         "branch": "0.001",
                         "branchn": "0.00004",
                         "cstore": "0.001",
                         "inorgn": "0.00004",
                         "metabsoil": "0.0",
                         "metabsoiln": "0.0",
                         "metabsurf": "0.0",
                         "metabsurfn": "0.0",
                         "nstore": "0.00004",
                         "passivesoil": "0.001",
                         "passivesoiln": "0.0004",
                         "prev_sma": "-999.9",
                         "root": "0.001",
                         "root_depth": "-9999.9",
                         "rootn": "0.00004",
                         "sapwood": "0.001",
                         "shoot": "0.001",
                         "shootn": "0.00004",
                         "slowsoil": "0.001",
                         "slowsoiln": "0.00004",
                         "stem": "0.001",
                         "stemn": "0.00004",
                         "stemnimm": "0.00004",
                         "stemnmob": "0.0",
                         "structsoil": "0.001",
                         "structsoiln": "0.00004",
                         "structsurf": "0.001",
                         "structsurfn": "0.00004",
                         "croot": "0.0",   # don't simulate coarse roots
                         "crootn": "0.0",  # don't simulate coarse roots

                         # parameters
                         "latitude": "35.9",
                         "intercep_frac": "0.15",
                         "max_intercep_lai": "3.0",
                         "albedo": "0.123",   # modis site avg
                         "finesoil": "0.51", # set based on silt+clay fractions of topsoil 0.42+0.09=0.5
                         "slamax": "4.4",     # Protocol [m2 kg-1 DW]
                         "sla": "4.4",        # Protocol [m2 kg-1 DW]
                         "slazero": "4.4",    # Protocol [m2 kg-1 DW]
                         "cfracts": "0.5",
                         "lai_closed": "0.5",  # I am effectively turning this feature off by setting it so low

                         #"c_alloc_fmax": "0.25",
                         #"c_alloc_fmin": "0.25",
                         #"c_alloc_rmax": "0.05",
                         #"c_alloc_rmin": "0.05",
                         #"c_alloc_bmax": "0.2",
                         #"c_alloc_bmin": "0.2",

                         #"c_alloc_fmax": "0.3",
                         #"c_alloc_fmin": "0.3",
                         #"c_alloc_rmax": "0.3",
                         #"c_alloc_rmin": "0.3",
                         #"c_alloc_bmax": "0.2",
                         #"c_alloc_bmin": "0.2",
                         #"c_alloc_cmax": "0.0", # turn off coarse roots!

                         "c_alloc_fmax": "0.35",
                         "c_alloc_fmin": "0.15",
                         "c_alloc_rmax": "0.35",
                         "c_alloc_rmin": "0.05",
                         "c_alloc_bmax": "0.1",
                         "c_alloc_bmin": "0.1",
                         "c_alloc_cmax": "0.0", # turn off coarse roots!


                         "fretrans": "0.5",
                         "rretrans": "0.0",
                         "bretrans": "0.0",
                         "wretrans": "0.0",
                         "ncwnewz": "0.003",
                         "ncwnew": "0.003",
                         "ncwimmz": "0.003",
                         "ncwimm": "0.003",
                         "ncbnewz": "0.003",
                         "ncbnew": "0.003",
                         "ncrfac": "0.8",
                         "ncmaxfyoung": "0.04",
                         "ncmaxfold": "0.04",
                         "ncmaxr": "0.03",
                         "retransmob": "0.0",
                         "fdecay": "0.59988",      # Protocol  [years-1]
                         "fdecaydry": "0.59988",   # Protocol
                         "rdecay": "0.33333",      # Protocol
                         "rdecaydry": "0.33333",   # Protocol
                         "bdecay": "0.02",         # No data, assuming 50 years
                         "wdecay": "0.02",
                         "crdecay": "0.00",  # turn off coarse roots!
                         "watdecaydry": "0.0",
                         "watdecaywet": "0.1",
                         "ligshoot": "0.24",      # Based on White et al. 2000 for ENF
                         "ligroot": "0.22",       # Based on White et al. 2000
                         "rateuptake": "2.2",           # set somewhat (very) arbitarly to get an LAI ~ 4.
                         "rateloss": "0.5",
                         "wcapac_root": "96.75",       # [mm] (FC (m3/m-3)-WP (m3/m-3)) * rooting_depth (mm) using derived values and depth from protocol, 750 mm (FC=0.164 - WP=0.035)
                         "wcapac_topsoil": "25.8",     # [mm] (FC (m3/m-3)-WP (m3/m-3)) * rooting_depth (mm) using derived values and depth from protocol, assuming 200 mm top soil following Corbeels 2005a (FC=0.164 - WP=0.035)


                         "ctheta_topsoil": "0.5",      # Derive based on soil type clay_loam
                         "ntheta_topsoil": "5.0",      # Derive based on soil type clay_loam
                         "ctheta_root": "0.4",         # Derive based on soil type clay
                         "ntheta_root": "3.0",         # Derive based on soil type clay
                         "topsoil_type": "clay_loam",
                         "rootsoil_type": "clay",
                         "measurement_temp": "25.0",
                         "dz0v_dh": "0.075", # However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                         "displace_ratio": "0.78",
                         "g1": "2.74",



                         #"jmaxna": "60.0",  # Original values Belinda had, I think based on Crous 2008, fig 2. Although those values I think are measured at 28 and 30 deg, the assumption being here that this is the same as 25 deg!
                         #"jmaxnb": "0.0",   # Original values Belinda had, I think based on Crous 2008, fig 2. Although those values I think are measured at 28 and 30 deg, the assumption being here that this is the same as 25 deg!
                         #"vcmaxna": "30.61",# Original values Belinda had, I think based on Crous 2008, fig 2. Although those values I think are measured at 28 and 30 deg, the assumption being here that this is the same as 25 deg!
                         #"vcmaxnb": "0.0",  # Original values Belinda had, I think based on Crous 2008, fig 2. Although those values I think are measured at 28 and 30 deg, the assumption being here that this is the same as 25 deg!
                         "vcmaxna": "22.29",
                         "vcmaxnb": "8.45",
                         "jv_slope": "1.86",
                         "jv_intercept": "0.0",



                         "sapturnover": "0.1",
                         "heighto": "4.826",
                         "htpower": "0.35",
                         "height0": "5.0",
                         "height1": "20.0",
                         "leafsap0": "8000.0",
                         "leafsap1": "3060.0",  # Duke protocol
                         "branch0": "5.61",
                         "branch1": "0.346",
                         "targ_sens": "0.5",
                         "density": "420.0",


                         # control
                         "adjust_rtslow": "false",  # priming, off
                         "alloc_model": "allometric",
                         "assim_model": "mate",
                         "calc_sw_params": "false",   #false=use fwp values, true=derive them
                         "deciduous_model": "false",
                         "disturbance": "false",
                         "exudation": "false",
                         "fixed_stem_nc": "true",
                         "fixleafnc": "false",
                         "grazing": "false",
                         "gs_model": "medlyn",
                         "model_optroot": "false",
                         "modeljm": "2",
                         "ncycle": "true",
                         "nuptake_model": "2",
                         "passiveconst": "false",
                         "print_options": "end",
                         "ps_pathway": "c3",
                         "respiration_model": "fixed",
                         "strfloat": "0",
                         "sw_stress_model": "1",  # Sands and Landsberg
                         "trans_model": "1",
                         "use_eff_nc": "0",
                         "use_leuning": "0",
                         "water_stress": "true",

                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY_SPIN + cfg_fname)

    if POST_INDUST == True:

        # run for 200 odd years post industrial with increasing co2/ndep
        # we are swapping forest params for grass params now
        # copy spunup base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, "%s_%s_model_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_spunup_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_spunup_adj" % (experiment_id, site)
        otag = "%s_%s_model_indust" % (experiment_id, site)
        mtag = "%s_met_data_industrial_to_present_1850_1983.csv" % (site)
        out_fn = itag + "_indust.out"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                         # git stuff
                         #"git_hash": str(git_revision),

                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),

                         # state - default C:N 25.
                         "age": "0.0",
                         "branch": "0.0",
                         "branchn": "0.0",
                         "canht": "0.79",  # Taken default C3grass value from JULES
                         "cstore": "0.001",
                         "nstore": "0.00004",
                         "croot": "0.0",   # don't simulate coarse roots
                         "crootn": "0.0",  # don't simulate coarse roots
                         "root": "0.001",
                         "rootn": "0.00004",
                         "sapwood": "0.0",
                         "shoot": "0.001",
                         "shootn": "0.00004",
                         "stem": "0.0",
                         "stemn": "0.0",
                         "stemnimm": "0.0",
                         "stemnmob": "0.0",
                         "nepsum": "0.0",
                         "nppsum": "0.0",

                         # parameters
                         "ligshoot": "0.09",  # Smith et al. 2000, GRASS
                         "ligroot": "0.22",   # Smith et al. 2000
                         "age": "1.0",
                         "slamax": "6.0",
                         "sla": "6.0",
                         "slazero": "6.0",
                         "cfracts": "0.5",
                         "lai_closed": "0.5",  # I am effectively turning this feature off by setting it so low
                         "c_alloc_fmax": "0.8",
                         "c_alloc_fmin": "0.2",
                         "c_alloc_rmax": "0.8",
                         "c_alloc_rmin": "0.2",
                         "c_alloc_bmax": "0.0",
                         "c_alloc_bmin": "0.0",
                         "c_alloc_cmax": "0.0", # turn off coarse roots!
                         "fretrans": "0.4",
                         "rretrans": "0.0",
                         "bretrans": "0.0",
                         "wretrans": "0.0",
                         "ncwnewz": "0.0",
                         "ncwnew": "0.0",
                         "ncwimmz": "0.0",
                         "ncwimm": "0.0",
                         "ncbnewz": "0.0",
                         "ncbnew": "0.0",
                         "ncrfac": "0.7",
                         "ncmaxfyoung": "0.035",
                         "ncmaxfold": "0.035",
                         "ncmaxr": "0.0287",
                         "retransmob": "0.0",
                         "fdecay": "1.0",
                         "fdecaydry": "1.0",
                         "rdecay": "1.0",
                         "rdecaydry": "1.0",
                         "bdecay": "0.0",
                         "wdecay": "0.0",
                         "watdecaydry": "0.0",
                         "watdecaywet": "0.1",
                         "crdecay": "0.00",  # turn off coarse roots!


                         "dz0v_dh": "0.10", # Taken default C3grass value from JULES
                         "displace_ratio": "0.64", #Jones 1992, pg. 67.
                         "z0h_z0m": "1.0",         # Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!!


                         "jmaxna": "62.0",  # assuming j = v * 2
                         "jmaxnb": "0.0",   # assuming no intercept
                         "vcmaxna": "31.0",   # C3 grasses - CLM4 tech doc, table 8.2, Oleson et al 2010, page 176
                         "vcmaxnb": "0.0",  # assuming no intercept

                         # control
                         "adjust_rtslow": "false",  # priming, off
                         "alloc_model": "grasses",
                         "assim_model": "mate",
                         "calc_sw_params": "false",   #false=use fwp values, true=derive them
                         "deciduous_model": "false",
                         "disturbance": "false",
                         "exudation": "false",
                         "fixed_stem_nc": "true",
                         "fixleafnc": "false",
                         "grazing": "false",
                         "gs_model": "medlyn",
                         "model_optroot": "false",
                         "modeljm": "1",
                         "nuptake_model": "2",
                         "ncycle": "true",
                         "passiveconst": "false",
                         "print_options": "end",
                         "ps_pathway": "c3",
                         "respiration_model": "fixed",
                         "strfloat": "0",
                         "sw_stress_model": "1",  # Sands and Landsberg
                         "trans_model": "1",
                         "use_eff_nc": "0",
                         "use_leuning": "0",
                         "water_stress": "true",

                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)


if __name__ == "__main__":

    experiment_id = "NCEAS"
    site = "DUKE"
    main(experiment_id, site, SPIN_UP=True, POST_INDUST=True)
