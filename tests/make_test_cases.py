#!/usr/bin/env python3
"""
Build the regression test cases from the Duke example.

The repository only ships daily forcing, so the daily Duke met is
disaggregated into a synthetic 30 min forcing (sinusoidal PAR over the day
length, a diurnal temperature/VPD cycle, rain spread evenly). This is a test
fixture, not realistic forcing. A drought variant has no rain after doy 90.

Usage: make_test_cases.py <example_dir> <out_dir> [nyears]
"""
import math
import os
import re
import sys

__author__ = "Martin De Kauwe"


def daylen(doy, ndays, lat):
    latr = math.radians(lat)
    sd = (-math.sin(math.radians(23.5)) *
          math.cos(2.0 * math.pi * (doy + 10.0) / ndays))
    a = math.sin(latr) * sd
    b = math.cos(latr) * math.cos(math.asin(sd))
    return 12.0 * (1.0 + (2.0 / math.pi) * math.asin(a / b))


def read_daily_met(fname, nyears):
    rows = [[float(x) for x in line.strip().split(",")]
            for line in open(fname) if not line.startswith("#")]
    first = int(rows[0][0])
    return [r for r in rows if int(r[0]) < first + nyears]


def write_30min_met(rows, fname, lat, drought=False):
    with open(fname, "w") as f:
        f.write("# synthetic 30 min Duke met (test fixture)\n")
        f.write("#year,doy,hod,rain,par,tair,tsoil,vpd,co2,ndep,nfix,wind,"
                "press\n")
        for r in rows:
            (yr, doy, tair, rain, tsoil, tam, tpm, tmin, tmax, tday, vpd_am,
             vpd_pm, co2, ndep, nfix, wind, press, wind_am, wind_pm, par_am,
             par_pm) = r
            yr, doy = int(yr), int(doy)
            ndays = 366 if yr % 4 == 0 else 365
            dl = daylen(doy, ndays, lat)
            sunrise, sunset = 12.0 - dl / 2.0, 12.0 + dl / 2.0
            if drought and doy >= 90:
                rain = 0.0

            # MJ PAR d-1 -> umol m-2 s-1 with a sinusoidal shape
            par_tot = (par_am + par_pm) * 1E6 * 4.57
            shape = [max(0.0, math.sin(math.pi * (h / 2.0 + 0.25 - sunrise) /
                                       dl))
                     if sunrise <= h / 2.0 + 0.25 <= sunset else 0.0
                     for h in range(48)]
            norm = sum(shape)
            for h in range(48):
                t = h / 2.0 + 0.25
                par = par_tot * shape[h] / norm / 1800.0 if norm > 0 else 0.0
                ta = tmin + (tmax - tmin) * 0.5 * \
                    (1.0 + math.sin(2.0 * math.pi * (t - 9.0) / 24.0))
                vpd = max(0.05, (vpd_am if t < 12 else vpd_pm) *
                          (0.3 + 0.7 * max(0.0, math.sin(math.pi * (t - 6.0) /
                                                         18.0))))
                f.write("%d,%d,%d,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6e,%.6e,"
                        "%.4f,%.4f\n" % (yr, doy, h, rain / 48.0, par, ta,
                                         tsoil, vpd, co2, ndep / 48.0,
                                         nfix / 48.0, wind, press))


def set_keys(txt, d):
    for key, val in d.items():
        pat = re.compile(r"^%s\s*=.*$" % key, re.M)
        if pat.search(txt):
            txt = pat.sub("%s = %s" % (key, val), txt)
        else:
            section = "[files]" if key.startswith("out_") else "[control]"
            if key == "root_psi_crit":
                section = "[params]"
            if key.startswith("soil_") and key not in ("soil_hydraulics",
                                                        "soil_drainage",
                                                        "soil_evap_model"):
                section = "[params]"
            txt = txt.replace(section, "%s\n%s = %s" % (section, key, val), 1)
    return txt


def main():
    example_dir, out_dir = sys.argv[1], sys.argv[2]
    nyears = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    os.makedirs(out_dir, exist_ok=True)

    base_cfg = os.path.join(example_dir, "params",
                            "NCEAS_DUKE_model_youngforest_amb.cfg")
    base = open(base_cfg).read()
    lat = float(re.search(r"^latitude\s*=\s*([-\d.]+)", base, re.M).group(1))
    rows = read_daily_met(os.path.join(example_dir, "met_data",
                                       "DUKE_met_data_amb_co2.csv"), nyears)
    write_30min_met(rows, os.path.join(out_dir, "met_30min.csv"), lat)
    write_30min_met(rows, os.path.join(out_dir, "met_30min_dry.csv"), lat,
                    drought=True)

    # daily example, as shipped
    daily = set_keys(base, {
        "met_fname": os.path.join(example_dir, "met_data",
                                  "DUKE_met_data_amb_co2.csv"),
        "out_fname": os.path.join(out_dir, "out_daily.csv"),
        "print_options": "daily"})
    open(os.path.join(out_dir, "daily.cfg"), "w").write(daily)

    # daily MATE with gs_opt (profit maximisation) on the bucket
    daily_gsopt = set_keys(base, {
        "met_fname": os.path.join(example_dir, "met_data",
                                  "DUKE_met_data_amb_co2.csv"),
        "out_fname": os.path.join(out_dir, "out_daily_gsopt.csv"),
        "print_options": "daily", "gs_model": "gs_opt",
        "nonstomatal_limitation": "false",
        "soil_hydraulics": "van_genuchten",
        "soil_b": "6.742", "soil_sathh": "0.22656875",
        "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"})
    open(os.path.join(out_dir, "daily_gsopt.cfg"), "w").write(daily_gsopt)

    # sub-daily cases
    cases = {
        "sd_bucket": {"water_balance": "bucket"},
        "sd_hyd": {"water_balance": "hydraulics"},
        "sd_hyd_dry": {"water_balance": "hydraulics", "met": "dry"},
        "sd_hyd_dry_cascade": {"water_balance": "hydraulics", "met": "dry",
                               "soil_drainage": "cascading"},
        "sd_hyd_dry_golden": {"water_balance": "hydraulics", "met": "dry",
                              "gs_opt_search": "golden"},
        "sd_hyd_dry_segmented": {"water_balance": "hydraulics", "met": "dry",
                                 "plant_segments": "segmented"},
        "sd_hyd_dry_jules_rz": {"water_balance": "hydraulics", "met": "dry",
                                "root_radial_resistance": "false",
                                "root_psi_crit": "-4.889"},
        # JULES FR-Pue soil parameters (b, sathh, satcon, sm_sat)
        "sd_hyd_dry_vg": {"water_balance": "hydraulics", "met": "dry",
                          "soil_hydraulics": "van_genuchten",
                          "bound_soil_psi": "true",
                          "dry_soil_correction": "true",
                          "soil_b": "6.742", "soil_sathh": "0.22656875",
                          "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"},
        "sd_hyd_vg_richards": {"water_balance": "hydraulics",
                               "soil_drainage": "richards",
                               "soil_hydraulics": "van_genuchten",
                               "soil_b": "6.742", "soil_sathh": "0.22656875",
                               "soil_satcon": "4.2228e-6",
                               "soil_sm_sat": "0.43648"},
        "sd_hyd_dry_vg_richards": {"water_balance": "hydraulics", "met": "dry",
                                   "soil_drainage": "richards",
                                   "soil_hydraulics": "van_genuchten",
                                   "bound_soil_psi": "true",
                                   "dry_soil_correction": "true",
                                   "soil_b": "6.742",
                                   "soil_sathh": "0.22656875",
                                   "soil_satcon": "4.2228e-6",
                                   "soil_sm_sat": "0.43648"},
        "sd_hyd_vg_richards_jules_esoil": {
            "water_balance": "hydraulics", "soil_drainage": "richards",
            "soil_hydraulics": "van_genuchten", "soil_evap_model": "jules",
            "soil_b": "6.742", "soil_sathh": "0.22656875",
            "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"},
        "sd_hyd_vg_richards_or_esoil": {
            "water_balance": "hydraulics", "soil_drainage": "richards",
            "soil_hydraulics": "van_genuchten", "soil_evap_model": "or",
            "soil_b": "6.742", "soil_sathh": "0.22656875",
            "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"},
        "sd_bucket_jules_esoil": {"water_balance": "bucket",
                                  "soil_evap_model": "jules"},
        "sd_hyd_vg_richards_jules_canopy": {
            "water_balance": "hydraulics", "soil_drainage": "richards",
            "soil_hydraulics": "van_genuchten", "soil_evap_model": "jules",
            "canopy_evap_model": "jules",
            "soil_b": "6.742", "soil_sathh": "0.22656875",
            "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"},
        "sd_bucket_jules_canopy": {"water_balance": "bucket",
                                   "canopy_evap_model": "jules"},
        "sd_hyd_dry_bc": {"water_balance": "hydraulics", "met": "dry",
                          "soil_hydraulics": "brooks_corey",
                          "soil_b": "6.742", "soil_sathh": "0.22656875",
                          "soil_satcon": "4.2228e-6", "soil_sm_sat": "0.43648"},
    }
    for name, opts in cases.items():
        met = "met_30min_dry.csv" if opts.pop("met", "") == "dry" else \
              "met_30min.csv"
        keys = {"met_fname": os.path.join(out_dir, met),
                "out_fname": os.path.join(out_dir, "out_%s.csv" % name),
                "out_subdaily_fname": os.path.join(out_dir,
                                                   "out_%s_30min.csv" % name),
                "sub_daily": "true", "print_options": "subdaily",
                "modeljm": "3", "vcmax": "60.0", "jmax": "100.0"}
        keys.update(opts)
        txt = set_keys(base, keys)
        open(os.path.join(out_dir, "%s.cfg" % name), "w").write(txt)

    print(" ".join(["daily", "daily_gsopt"] + list(cases)))


if __name__ == "__main__":
    main()
