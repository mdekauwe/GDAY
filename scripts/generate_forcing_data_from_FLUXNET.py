#!/usr/bin/env python

"""
Create a G'DAY met forcing file from a FLUXNET file.

FLUXNET files first need to be created via Anna Ukkola's package.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (12.09.2018)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import csv
import math
import numpy as np
from datetime import date
import calendar
import pandas as pd
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import datetime as dt
import netCDF4 as nc
import xarray as xr

class CreateMetData(object):

    def __init__(self, fdir, site, daily=True, num_yrs=20, co2_spinup=285.):

        self.fdir = fdir
        self.site = site
        self.fname = os.path.join(self.fdir, "%sOzFlux2.0_met.nc" % (site))
        self.spinup_ofname = "%s_met_spinup.csv" % (site)
        self.forcing_ofname = "%s_met_forcing.csv" % (site)
        self.daily = daily

        if self.daily:
            self.ovar_names = ['#year', 'doy', 'tair', 'rain', 'tsoil',
                               'tam', 'tpm', 'tmin', 'tmax', 'tday', 'vpd_am',
                               'vpd_pm', 'co2', 'ndep', 'nfix', 'wind', 'press',
                               'wind_am', 'wind_pm', 'par_am', 'par_pm']
            self.ounits = ['#--', '--', 'degC', 'mm/d', 'degC','degC', 'degC',
                           'degC','degC', 'degC', 'kPa', 'kPa', 'ppm', 't/ha/d',
                           't/ha/d', 'm/s', 'kPa', 'm/s', 'm/s', 'mj/m2/d',
                           'mj/m2/d']
        else:
            self.ovar_names = ['#year', 'doy', 'hod', 'rain', 'par', 'tair',
                               'tsoil', 'vpd', 'co2', 'ndep', 'nfix', 'wind',
                               'press']
            self.ounits = ['#--', '--', '--', 'mm/30min', 'umol/m2/s','degC',
                           'degC', 'kPa', 'ppm', 't/ha/30min', 'm/s','kPa']
        self.lat = -999.9
        self.lon = -999.9
        self.num_yrs = num_yrs
        self.co2_spinup = co2_spinup

        # unit conversions
        self.SEC_TO_HLFHR = 1800.0 # sec/half hr
        self.SEC_TO_HR = 3600.0 # sec/ hr
        self.PA_TO_KPA = 0.001
        self.J_TO_MJ = 1.0E-6
        self.K_to_C = 273.15
        self.G_M2_to_TONNES_HA = 0.01
        self.MM_S_TO_MM_HR = 3600.0
        self.SW_2_PAR = 2.3
        self.J_TO_UMOL = 4.57
        self.UMOL_TO_J = 1.0 / self.J_TO_UMOL

    def main(self, vary_ndep=False, ndep=-999.9):

        df = self.read_nc_file()

        start_yr = df.index.year[0]
        end_yr = df.index.year[-1]
        yr_sequence = self.get_year_sequence(start_yr, end_yr)

        self.write_met_file(df, spinup=True, yr_sequence=yr_sequence,
                            vary_co2=False, co2_data=self.co2_spinup,
                            vary_ndep=False, ndep_data=-999.9, vary_nfix=False,
                            nfix_data=-999.9)
        self.write_met_file(df, spinup=False, yr_sequence=None, vary_co2=True,
                            co2_data=None, vary_ndep=False, ndep_data=-999.9,
                            vary_nfix=False, nfix_data=-999.9)

    def write_met_file(self, df, spinup=True, yr_sequence=None, vary_co2=False,
                       co2_data=None, vary_ndep=False, ndep_data=None,
                       vary_nfix=False, nfix_data=None):

        if spinup == False:
            yr_sequence = np.unique(df.year)

        start_sim = yr_sequence[0]
        end_sim = yr_sequence[-1]
        year = str(start_sim)

        try:
            if spinup:
                if os.path.isfile(self.spinup_ofname):
                    os.remove(self.spinup_ofname)
                ofp = open(self.spinup_ofname, 'w')
            else:
                if os.path.isfile(self.forcing_ofname):
                    os.remove(self.forcing_ofname)
                ofp = open(self.forcing_ofname, 'w')
            wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                            escapechar=None, dialect='excel')
            if self.daily:
                wr.writerow(['# %s daily met forcing' % (site)])
            else:
                wr.writerow(['# %s sub-daily met forcing' % (site)])
            wr.writerow(['# Data from %s-%s' % (start_sim, end_sim)])
            wr.writerow(['# Created by Martin De Kauwe: %s' % date.today()])
            wr.writerow([var for i, var in enumerate(self.ounits)])
            wr.writerow([var for i, var in enumerate(self.ovar_names)])

        except IOError:
            raise IOError('Could not write met file: %s' % \
                          self.spinup_ofname)

        # Account for hourly vs half-hourly files.
        diff = df.index.minute[1] - df.index.minute[0]

        for i, yr in enumerate(yr_sequence):
            days = np.unique(df[df.year == yr].doy)
            for j, doy in enumerate(days):
                days_data = df[(df.year == yr) & (df.doy == doy)]

                # otherwise we can't index 0-47 for HOD
                days_data = days_data.reset_index()

                if self.daily:
                    self.daily_unpacking(yr, doy, days_data, wr, diff,
                                         vary_co2, vary_ndep, vary_nfix,
                                         co2_data, ndep_data, nfix_data)
                else:
                    self.sub_daily_unpacking(yr, doy, days_data, wr, diff,
                                             vary_co2, vary_ndep, vary_nfix,
                                             co2_data, ndep_data, nfix_data)
        ofp.close()

    def read_nc_file(self):
        """ Build a dataframe from the netcdf outputs """

        ds = xr.open_dataset(self.fname)


        self.lat = ds.latitude.values[0][0]
        self.lon = ds.longitude.values[0][0]

        # W/m2, deg K, mm/s, kg/kg, m/s, Pa, ppmw
        vars_to_keep = ['SWdown','Tair','Rainf','Qair','Wind','PSurf','CO2air']
        df = ds[vars_to_keep].squeeze(dim=["x","y","z"],
                                      drop=True).to_dataframe()

        # PALS-style netcdf is missing the first (half)hour timestamp and has
        # one extra from the next year, i.e. everything is shifted, so we need
        # to fix this. We will duplicate the
        # first hour interval and remove the last
        time_idx = df.index
        diff = df.index.minute[1] - df.index.minute[0]
        if diff == 0:
            time_idx = time_idx.shift(-1, freq='H')
            df = df.shift(-1, freq='H')
        else:
            time_idx = time_idx.shift(-1, freq='30min')
            df = df.shift(-1, freq='30min')

        df = df.reindex(time_idx)
        df['year'] = df.index.year
        df['doy'] = df.index.dayofyear

        df["PAR"] = df.SWdown * self.SW_2_PAR

        return df

    def get_year_sequence(self, start_yr, end_yr):

        yrs = np.arange(start_yr, end_yr+1)
        sequence = []
        i = 0
        while len(sequence) != self.num_yrs:
            sequence.append(yrs[i])
            i += 1
            if i == len(yrs):
                i = 0

        return sequence

    def qair_to_vpd(self, qair, tair, press):

        # convert back to Pa
        press /= self.PA_TO_KPA

        # saturation vapor pressure
        es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

        # vapor pressure
        ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

        vpd = (es - ea) * self.PA_TO_KPA

        return vpd

    def sub_daily_unpacking(self, yr, doy, days_data,  wr, diff, vary_co2,
                            vary_ndep, vary_nfix, co2_data, ndep_data,
                            nfix_data):

        for hod in range(len(days_data)):

            # mm/sec -> mm/30 min
            # hour gap i.e. Tumba
            if diff == 0:
                # need to split the rain in two & spread over 2 half hrs
                rain = days_data.Rainf[hod] * self.MM_S_TO_MM_HR / 2.0
            else:
                rain = days_data.Rainf[hod] * self.MM_S_TO_MM_HR

            par = days_data.SWdown[hod] * self.SW_2_PAR
            if par < 0.0:
                par = 0.0
            tair = days_data.Tair[hod] - self.K_to_C
            tsoil = np.mean(days_data.Tair) - self.K_to_C
            qair = days_data.Qair[hod]

            # co2 -> [ppm] Daily mean value
            if vary_co2:
                co2 = days_data.CO2air[hod]
            else:
                co2 = co2_data

            # g/m2/30 min -> t/ha/30min
            if vary_ndep:
                ndep = ndep_data[hod] * self.G_M2_to_TONNES_HA
            else:
                ndep = ndep_data

            if vary_nfix:
                nfix = nfix_data[hod] * self.G_M2_to_TONNES_HA
            else:
                nfix = nfix_data

            wind = days_data.Wind[hod]
            if wind <= 0.0:
                wind = 0.1 # set v.small speed but not zero
            press = days_data.PSurf[hod] * self.PA_TO_KPA
            vpd = self.qair_to_vpd(qair, tair, press)
            vpd = np.where(vpd < 0.05, 0.05, vpd)

            if diff == 0:
                # Need to write two rows as the data is hourly and we
                # want 30 min
                wr.writerow([yr, doy, hod*2., rain, par, tair, tsoil, \
                             vpd, co2, ndep, nfix, wind, press])
                wr.writerow([yr, doy, hod*2.+0.5, rain, par, tair, \
                             tsoil, vpd, co2, ndep, nfix, wind, press])
            else:
                wr.writerow([yr, doy, hod, rain, par, tair, tsoil, \
                             vpd, co2, ndep, nfix, wind, press])

    def daily_unpacking(self, yr, doy, days_data,  wr, diff, vary_co2,
                        vary_ndep, vary_nfix, co2_data, ndep_data, nfix_data):

        if diff == 0:
            morning = days_data.iloc[0:12]
            morning = morning[morning.PAR >= 5.0]
            afternoon = days_data.iloc[12:24]
            afternoon = afternoon[afternoon.PAR >= 5.0]
        else:
            morning = days_data.iloc[0:24]
            morning = morning[morning.PAR >= 5.0]
            afternoon = days_data.iloc[24:48]
            afternoon = afternoon[afternoon.PAR >= 5.0]
        day_light = days_data[days_data.PAR >= 5.0]

        tair = np.mean(day_light.Tair - self.K_to_C)
        tam = np.mean(morning.Tair - self.K_to_C)
        tpm = np.mean(afternoon.Tair - self.K_to_C)
        tsoil = np.mean(days_data.Tair) - self.K_to_C

        # daytime min/max temp
        tmin = np.min(days_data.Tair - self.K_to_C)
        tmax = np.max(days_data.Tair - self.K_to_C)
        tday = np.mean(days_data.Tair - self.K_to_C)

        # set v.small speed but not zero
        wind = np.mean(day_light.Wind)
        wind = np.where(wind <= 0.0, 0.1, wind)
        wind_am = np.mean(morning.Wind)
        wind_am = np.where(wind_am <= 0.0, 0.1, wind_am)
        wind_pm = np.mean(afternoon.Wind)
        wind_pm = np.where(wind_pm <= 0.0, 0.1, wind_pm)

        press = np.mean(day_light.PSurf) * self.PA_TO_KPA
        press_am = np.mean(morning.PSurf) * self.PA_TO_KPA
        press_pm = np.mean(afternoon.PSurf) * self.PA_TO_KPA
        qair_am = np.mean(morning.Qair)
        qair_pm = np.mean(afternoon.Qair)

        vpd_am = self.qair_to_vpd(qair_am, tam, press_am)
        vpd_am = np.where(vpd_am < 0.05, 0.05, vpd_am)
        vpd_pm = self.qair_to_vpd(qair_pm, tpm, press_pm)
        vpd_pm = np.where(vpd_am < 0.05, 0.05, vpd_pm)

        # convert PAR [umol m-2 s-1] -> mj m-2 30min-1
        if diff == 0:
            conv = self.UMOL_TO_J * self.J_TO_MJ * self.SEC_TO_HR
        else:
            conv = self.UMOL_TO_J * self.J_TO_MJ * self.SEC_TO_HLFHR
        par_am = np.sum(morning.PAR * conv)
        par_pm = np.sum(afternoon.PAR * conv)

        # rain -> mm/30 min
        # Don't need the conversion!
        # coversion 1800 seconds to half hours and summed gives day value.
        # Going to use the whole day including the night data
        #rain = max(0.0, np.sum(days_data.Precip * 1800.))
        rain = max(0.0, np.sum(days_data.Rainf))

        # co2 -> [ppm] Daily mean value
        if vary_co2:
            co2 = np.mean(day_light.CO2air)
        else:
            co2 = co2_data

        # g/m2/30 min -> t/ha/30min
        if vary_ndep:
            ndep = ndep_data * self.G_M2_to_TONNES_HA
        else:
            ndep = ndep_data

        if vary_nfix:
            nfix = nfix_data * self.G_M2_to_TONNES_HA
        else:
            nfix = nfix_data

        wr.writerow([yr, doy, tair, rain, tsoil, tam, tpm, tmin, tmax,\
                             tday, vpd_am, vpd_pm, co2, ndep, nfix, wind, \
                             press, wind_am, wind_pm, par_am, par_pm])


if __name__ == "__main__":

    fdir = "/Users/%s/research/OzFlux" % (os.getlogin())
    site = "Tumbarumba"

    # Daily file
    C = CreateMetData(fdir, site, daily=True)
    C.main()
