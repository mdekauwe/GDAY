#!/usr/bin/env python

""" 
Series of unit-tests for GDAY code to pass at installation time, or
when new code has been changed. 
"""

import os
import sys
import numpy as np
import unittest
from math import exp, sqrt, sin, pi
from gday.mate import MateC3, MateC4
from gday.file_parser import read_met_forcing
import gday.default_control as control
import gday.default_files as files
import gday.default_params as params
import gday.default_fluxes as fluxes
import gday.default_state as state
from gday.water_balance import WaterBalance, SoilMoisture

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.012.2014)"
__email__   = "mdekauwe@gmail.com"



def setup_metdata(day):
    #met_fname = "../example/met_data/DUKE_met_data_amb_co2.csv"
    #met_data = read_met_forcing(met_fname, met_header=4)
    #print met_data['tam'][day] 
    #print met_data['tpm'][day] 
    #print met_data['tair'][day] 
    #print met_data['sw_rad'][day]
    #print met_data['sw_rad_am'][day] 
    #print met_data['sw_rad_pm'][day] 
    #print met_data["rain"][day] 
    #print met_data['vpd_am'][day] 
    #print met_data['vpd_pm'][day]
    #print met_data['vpd_avg'][day]
    #print met_data['wind_am'][day]
    #print met_data['wind_pm'][day]
    #print met_data['wind'][day]
    #print met_data["co2"][day]
    #print met_data['atmos_press'][day] 

    met_data = {}
    met_data['tam'] = {} 
    met_data['tam'][day] = {} 
    met_data['tpm'] = {} 
    met_data['tpm'][day] = {} 
    met_data['tair'] = {} 
    met_data['tair'][day] = {} 
    
    met_data['sw_rad'] = {} 
    met_data['sw_rad'][day] = {} 
    met_data['sw_rad_am'] = {} 
    met_data['sw_rad_am'][day] = {} 
    met_data['sw_rad_pm'] = {} 
    met_data['sw_rad_pm'][day] = {} 
    
    met_data['rain'] = {} 
    met_data['rain'][day] = {} 
    
    met_data['vpd_avg'] = {} 
    met_data['vpd_avg'][day] = {} 
    met_data['vpd_am'] = {} 
    met_data['vpd_am'][day] = {} 
    met_data['vpd_pm'] = {} 
    met_data['vpd_pm'][day] = {} 
    
    met_data['co2'] = {} 
    met_data['co2'][day] = {} 
    
    met_data['atmos_press'] = {} 
    met_data['atmos_press'][day] = {} 
    
    met_data['wind_am'] = {} 
    met_data['wind_am'][day] = {} 
    met_data['wind_pm'] = {} 
    met_data['wind_pm'][day] = {} 
    met_data['wind'] = {} 
    met_data['wind'][day] = {} 
    met_data['par'] = {} 
    met_data['par'][day] = {} 
    
    # Add values
    met_data['tam'][day] = 10.0
    met_data['tpm'][day] = 25.0
    met_data['tair'][day] = (10.0 + 25.0) / 2.0
    met_data['sw_rad'][day] = 20.5696188
    met_data['sw_rad_am'][day] = 11.9967552
    met_data['sw_rad_pm'][day] = 8.5728636
    met_data["rain"][day] = 9.0
    met_data['vpd_am'][day] = 0.3
    met_data['vpd_pm'][day] = 0.8
    met_data['vpd_avg'][day] = 0.6
    met_data['wind_am'][day] = 3.65
    met_data['wind_pm'][day] = 6.1
    met_data['wind'][day] = 4.8
    met_data["co2"][day] = 380.0
    met_data['atmos_press'][day] = 99.7
    met_data['par'][day] = 4674600.0
    
    return met_data
    
def setup_params():

    params.a1 = 0.0
    params.ac = 0.5
    params.actncmax = 0.333333
    params.actncmin = 0.066667
    params.adapt = 0.012
    params.ageold = 10000.0
    params.ageyoung = 0.0
    params.albedo = 0.123
    params.alpha_c4 = 0.06
    params.alpha_j = 0.26
    params.b_root = None
    params.b_topsoil = None
    params.bdecay = 0.02
    params.branch0 = 5.61
    params.branch1 = 0.346
    params.bretrans = 0.0
    params.burn_specific_yr = None
    params.c_alloc_bmax = 0.2
    params.c_alloc_bmin = 0.2
    params.c_alloc_cmax = 0.0
    params.c_alloc_fmax = 0.3
    params.c_alloc_fmin = 0.3
    params.c_alloc_rmax = 0.3
    params.c_alloc_rmin = 0.3
    params.cfracts = 0.5
    params.crdecay = 0.0
    params.cretrans = 0.0
    params.croot0 = 0.34
    params.croot1 = 0.84
    params.ctheta_root = 0.4
    params.ctheta_topsoil = 0.5
    params.cue = 0.5
    params.d0 = 0.0
    params.d0x = 0.35
    params.d1 = 0.0
    params.delsj = 644.4338
    params.density = 420.0
    params.direct_frac = 0.5
    params.displace_ratio = 0.78
    params.disturbance_doy = 1.0
    params.dz0v_dh = 0.075
    params.eac = 79430.0
    params.eag = 37830.0
    params.eaj = 43790.0
    params.eao = 36380.0
    params.eav = 51560.0
    params.edj = 200000.0
    params.faecescn = 25.0
    params.faecesn = 0.0
    params.fdecay = 0.59988
    params.fdecaydry = 0.59988
    params.fhw = 0.8
    params.finesoil = 0.51
    params.fracfaeces = 0.3
    params.fracteaten = 0.5
    params.fractosoil = 0.85
    params.fractup_soil = 0.5
    params.fretrans = 0.5
    params.g1 = 2.74
    params.gamstar25 = 42.75
    params.growth_efficiency = 0.7
    params.height0 = 5.0
    params.height1 = 25.0
    params.heighto = 4.826
    params.htpower = 0.35
    params.hurricane_doy = None
    params.hurricane_yr = None
    params.intercep_frac = 0.15
    params.jmax = -999.9
    params.jmaxna = 41.4594
    params.jmaxnb = 0.0
    params.kc25 = 404.9
    params.kdec1 = 3.965571
    params.kdec2 = 14.61
    params.kdec3 = 4.904786
    params.kdec4 = 18.262499
    params.kdec5 = 7.305
    params.kdec6 = 0.198279
    params.kdec7 = 0.006783
    params.kext = 0.5
    params.knl = 0.01
    params.ko25 = 278400.0
    params.kq10 = 0.08
    params.kr = 0.5
    params.lai_closed = 0.5
    params.latitude = 35.9
    params.leafsap0 = 10000.0
    params.leafsap1 = 3060.0
    params.ligfaeces = 0.25
    params.ligroot = 0.22
    params.ligshoot = 0.24
    params.liteffnc = 0.0
    params.max_intercep_lai = 3.0
    params.measurement_temp = 25.0
    params.ncbnew = 0.003
    params.ncbnewz = 0.003
    params.nccnew = 0.003
    params.nccnewz = 0.003
    params.ncmaxfold = 0.06
    params.ncmaxfyoung = 0.06
    params.ncmaxr = 0.03
    params.ncrfac = 0.8
    params.ncwimm = 0.003
    params.ncwimmz = 0.003
    params.ncwnew = 0.003
    params.ncwnewz = 0.003
    params.nf_crit = 0.015
    params.nf_min = 0.005
    params.nmax = 0.24
    params.nmin = 0.95
    params.nmin0 = 0.0
    params.nmincrit = 2.0
    params.ntheta_root = 3.0
    params.ntheta_topsoil = 5.0
    params.nuptakez = 0.0
    params.oi = 205000.0
    params.passivesoilnz = 1.0
    params.passivesoilz = 1.0
    params.passncmax = 0.142857
    params.passncmin = 0.1
    params.prescribed_leaf_NC = 0.03
    params.previous_ncd = 35.0
    params.prime_y = 0.0
    params.prime_z = 0.0
    params.psi_sat_root = None
    params.psi_sat_topsoil = None
    params.qs = 1.0
    params.r0 = 0.1325
    params.rateloss = 0.5
    params.rateuptake = 2.0
    params.rdecay = 0.33333
    params.rdecaydry = 0.33333
    params.retransmob = 0.0
    params.return_interval = 10.0
    params.rfmult = 1.0
    params.root_exu_CUE = None
    params.rooting_depth = 750.0
    params.rootsoil_type = "clay"
    params.rretrans = 0.0
    params.sapturnover = 0.1
    params.sla = 4.4
    params.slamax = 4.4
    params.slazero = 4.4
    params.slowncmax = 0.066666
    params.slowncmin = 0.025
    params.store_transfer_len = None
    params.structcn = 150.0
    params.structrat = 0.0
    params.targ_sens = 0.5
    params.theta = 0.7
    params.theta_sat_root = None
    params.theta_sat_topsoil = None
    params.topsoil_depth = 350.0
    params.topsoil_type = "clay_loam"
    params.vcmax = -999.9
    params.vcmaxna = 22.29
    params.vcmaxnb = 8.45
    params.watdecaydry = 0.0
    params.watdecaywet = 0.1
    params.wcapac_root = 96.75
    params.wcapac_topsoil = 25.8
    params.wdecay = 0.02
    params.wetloss = 0.5
    params.wretrans = 0.0
    params.z0h_z0m = 1.0 
    
    return params

def testPhotosynthesis(ps_pathway=None, debug=False):
    
    # Setup stuff
    day = 100
    params = setup_params()
    met_data = setup_metdata(day)
    
    state.wtfac_root = 0.8
    daylen = 12.0
    state.ncontent = 4.5
    state.lai = 3.0
    state.fipar = (1.0 - exp(-params.kext * state.lai))
    
    if ps_pathway == "C3":
        M = MateC3(control, params, state, fluxes, met_data)
        M.calculate_photosynthesis(day, daylen)
    elif ps_pathway == "C4":
        M = MateC4(control, params, state, fluxes, met_data)
        M.calculate_photosynthesis(day, daylen)
    if debug:
        print "Mate"
        print fluxes.gpp_gCm2
        print fluxes.npp_gCm2
        print fluxes.gpp_am
        print fluxes.gpp_pm
   
     
def testWater(ps_pathway=None, debug=False):
    # Setup stuff
    day = 100
    params = setup_params()
    met_data = setup_metdata(day)
    testPhotosynthesis(ps_pathway="C3")
    
    project_day = 100
    daylen = 12.0
    state.ncontent = 4.5
    state.lai = 3.0
    state.fipar = (1.0 - exp(-params.kext * state.lai))
    state.pawater_root = 50.2620995991
    state.pawater_topsoil = 10.734719666
   
    wb = WaterBalance(control, params, state, fluxes, met_data)
    sm = SoilMoisture(control, params, state, fluxes)
    (state.wtfac_topsoil, wtfac_root) = sm.calculate_soil_water_fac()
    wb.calculate_water_balance(project_day, daylen)
    
    if debug:
        print "Water"   
        print state.wtfac_topsoil
        print state.wtfac_root
        print fluxes.et
        print fluxes.soil_evap
        print fluxes.transpiration
        print fluxes.erain
        print fluxes.interception
        print fluxes.runoff
        print fluxes.gs_mol_m2_sec
        print fluxes.ga_mol_m2_sec
    
class GdayTests(unittest.TestCase):
    
    print 
    print "---------------------"
    print
    
    def testMateC3(self):
        print "Testing Photosynthesis - C3"
        print 
        testPhotosynthesis(ps_pathway="C3")
        def test_total_gpp(self):
            # Values pre-calculated
            gpp_gCm2 = 1.85043009942
            self.assertAlmostEqual(gpp_gCm2, fluxes.gpp_gCm2)
    
        def test_total_npp(self):
            # Values pre-calculated
            npp_gCm2 = 0.925215049712
            self.assertAlmostEqual(npp_gCm2, fluxes.npp_gCm2)
    
        def test_gpp_am(self):    
            # Values pre-calculated
            gpp_am = 1.00246350983
            self.assertAlmostEqual(gpp_am, fluxes.gpp_am)
            
    
        def test_gpp_pm(self):       
            # Values pre-calculated
            correct_value = 0.847966589596
            self.assertAlmostEqual(correct_value, fluxes.gpp_pm)
      
    def testMateC4(self):
        print "Testing Photosynthesis - C4"
        print 
        testPhotosynthesis(ps_pathway="C4")
        def test_total_gpp(self):
            # Values pre-calculated
            correct_value = 1.86246931786
            self.assertAlmostEqual(correct_value, fluxes.gpp_gCm2)
    
        def test_total_npp(self):
            # Values pre-calculated
            correct_value = 0.931234658932
            self.assertAlmostEqual(correct_value, fluxes.npp_gCm2)
    
        def test_gpp_am(self):    
            # Values pre-calculated
            correct_value = 0.959254213558
            self.assertAlmostEqual(correct_value, fluxes.gpp_am)
           
    
        def test_gpp_pm(self):       
            # Values pre-calculated
            correct_value = 0.903215104306     
            self.assertAlmostEqual(correct_value, fluxes.gpp_pm)
    
    def testWaterBalance3(self):
        print "Testing Water Balance"
        print 
        testWater(ps_pathway="C3")   
        def test_wtfac_topsoil(self):
            # Values pre-calculated
            correct_value = 0.315219884377
            self.assertAlmostEqual(correct_value, state.wtfac_topsoil)
        def test_wtfac_root(self):
            # Values pre-calculated
            correct_value = 0.8
            self.assertAlmostEqual(correct_value, state.wtfac_root)
        def test_et(self):
            # Values pre-calculated
            correct_value = 1.96794847661
            self.assertAlmostEqual(correct_value, fluxes.et)
        def test_soil_evap(self):
            # Values pre-calculated
            correct_value = 0.34793507562
            self.assertAlmostEqual(correct_value, fluxes.soil_evap)
        def test_transpiration(self):
            # Values pre-calculated
            correct_value = 0.270013400995
            self.assertAlmostEqual(correct_value, fluxes.transpiration)
        def test_erain(self):
            # Values pre-calculated
            correct_value = 7.65
            self.assertAlmostEqual(correct_value, fluxes.erain)
        def test_interception(self):
            # Values pre-calculated
            correct_value = 1.35
            self.assertAlmostEqual(correct_value, fluxes.interception)
        def test_runoff(self):
            # Values pre-calculated
            correct_value = 0.0
            self.assertAlmostEqual(correct_value, fluxes.runoff)
        def test_ga(self):
            # Values pre-calculated
            correct_value = 0.0644936068115
            self.assertAlmostEqual(correct_value, fluxes.ga_mol_m2_sec)
        def test_gs(self):
            # Values pre-calculated
            correct_value = 29.0257955939
            self.assertAlmostEqual(correct_value, fluxes.gs_mol_m2_sec)
   

    
    
        
if __name__ == "__main__":
    
    
    #testPhotosynthesis(ps_pathway="C3", debug=True)
    #testWater(ps_pathway="C3", debug=True)
    
    unittest.main()