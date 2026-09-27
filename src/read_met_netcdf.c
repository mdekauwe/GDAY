/* ============================================================================
* Read met forcing from a netCDF file in the PLUMBER2 / ALMA convention (the
* same files JULES reads), e.g. FR-Pue_2000-2014_FLUXNET2015_Met.nc
*
* Variables used (time, y, x) with y = x = 1:
*   Tair (K), SWdown (W m-2), Precip (kg m-2 s-1), Qair (kg kg-1),
*   Psurf (Pa), Wind (m s-1), CO2air (ppm, optional),
*   LWdown (W m-2, optional, sub-daily only)
*
* The sub-daily model uses the timesteps directly (30 min data required).
* For the daily model the timesteps are aggregated into the daily am/pm
* forcing in the same way as scripts/generate_forcing_data_from_FLUXNET.py:
* daylight (PAR >= 5 umol m-2 s-1) timesteps before/after noon.
*
* Not in PLUMBER2 files: soil temperature (daily mean air temperature is
* used), N deposition & fixation (params nc_ndep, nc_nfix, t N ha-1 yr-1).
*
* Prescribed LAI (control prescribed_lai) is read from lai_fname if set
* (e.g. the JULES MODIS file, variable lai_var, dimension lai_pft_index),
* otherwise from the met file's LAI variable. The LAI file either has the
* met file's timesteps, or is daily ("days/seconds since" time units), in
* which case each day's value is held over that day's timesteps.
*
* Only built with HAVE_NETCDF (see the Makefile).
* =========================================================================== */
#include "read_met_file.h"

#ifdef HAVE_NETCDF
#include <netcdf.h>

#define NC_CHECK(stat, what) { int _s = (stat); if (_s != NC_NOERR) { \
    fprintf(stderr, "netCDF error (%s): %s\n", what, nc_strerror(_s)); \
    exit(EXIT_FAILURE); } }

static long days_from_civil(int y, int m, int d) {
    /* days since 1970-01-01 (proleptic Gregorian), H. Hinnant's algorithm */
    y -= m <= 2;
    long era = (y >= 0 ? y : y - 399) / 400;
    long yoe = y - era * 400;
    long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return (era * 146097 + doe - 719468);
}

static void civil_from_days(long z, int *y, int *doy) {
    /* year and day of year (1-based) from days since 1970-01-01 */
    int  yy, m, d;
    z += 719468;
    long era = (z >= 0 ? z : z - 146096) / 146097;
    long doe = z - era * 146097;
    long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    yy = (int)(yoe + era * 400);
    long dy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    long mp = (5 * dy + 2) / 153;
    d = (int)(dy - (153 * mp + 2) / 5 + 1);
    m = (int)(mp < 10 ? mp + 3 : mp - 9);
    yy += (m <= 2);
    *y = yy;
    *doy = (int)(days_from_civil(yy, m, d) - days_from_civil(yy, 1, 1)) + 1;
}

static double *alloc_array(size_t, const char *);

static double *read_nc_var(int ncid, const char *name, size_t n,
                           int required) {
    /* read a (time[, y, x]) variable as double; NULL if absent & optional */
    int    varid, stat;
    size_t i;
    float  fill = -9999.0;
    float *buf;
    double *out;

    stat = nc_inq_varid(ncid, name, &varid);
    if (stat != NC_NOERR) {
        if (required) {
            fprintf(stderr, "netCDF met file has no variable %s\n", name);
            exit(EXIT_FAILURE);
        }
        return (NULL);
    }
    if ((buf = malloc(n * sizeof(float))) == NULL ||
        (out = malloc(n * sizeof(double))) == NULL) {
        fprintf(stderr, "Error allocating space for %s\n", name);
        exit(EXIT_FAILURE);
    }
    NC_CHECK(nc_get_var_float(ncid, varid, buf), name);
    nc_get_att_float(ncid, varid, "_FillValue", &fill);
    for (i = 0; i < n; i++) {
        if (buf[i] == fill || isnan(buf[i])) {
            fprintf(stderr, "netCDF met variable %s is missing at timestep "
                    "%zu, the forcing must be gap filled\n", name, i);
            exit(EXIT_FAILURE);
        }
        out[i] = (double)buf[i];
    }
    free(buf);

    return (out);
}

static long *read_lai_days(int ncid, size_t len) {
    /* day (since 1970-01-01) of each record of a daily LAI file */
    int     varid, y, mo, d, h = 0, mi = 0, s = 0;
    size_t  i;
    double *t, scale;
    long   *day;
    char    units[NC_MAX_NAME + 1], what[16];

    NC_CHECK(nc_inq_varid(ncid, "time", &varid), "LAI time variable");
    memset(units, 0, sizeof(units));
    NC_CHECK(nc_get_att_text(ncid, varid, "units", units), "LAI time units");
    if (sscanf(units, "%15s since %d-%d-%d %d:%d:%d", what, &y, &mo, &d, &h,
               &mi, &s) < 4) {
        fprintf(stderr, "Can't parse LAI time units '%s'\n", units);
        exit(EXIT_FAILURE);
    }
    if (strcmp(what, "days") == 0) {
        scale = 1.0;
    } else if (strcmp(what, "seconds") == 0) {
        scale = 1.0 / 86400.0;
    } else {
        fprintf(stderr, "LAI time units must be days or seconds since ...\n");
        exit(EXIT_FAILURE);
    }
    t = alloc_array(len, "LAI time");
    day = malloc(len * sizeof(long));
    NC_CHECK(nc_get_var_double(ncid, varid, t), "LAI time");
    for (i = 0; i < len; i++) {
        day[i] = days_from_civil(y, mo, d) +
                 (long)floor((h * 3600.0 + mi * 60.0 + s) / 86400.0 +
                             t[i] * scale + 1E-9);
        if (i > 0 && day[i] != day[i-1] + 1) {
            fprintf(stderr, "Daily LAI file must have consecutive days "
                    "(record %zu)\n", i);
            exit(EXIT_FAILURE);
        }
    }
    free(t);

    return (day);
}

static double *read_nc_lai(control *c, size_t n, int met_ncid,
                           const long *met_day) {
    /* prescribed LAI, from lai_fname (pft dimension) or the met file.
       met_day: day (since 1970-01-01) of each met timestep */
    int     ncid, varid, ndims, dimids[NC_MAX_VAR_DIMS], d, idim_time = -1;
    size_t  start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS], len, i, nlai = n;
    long    k, *lai_day = NULL;
    double *out, *daily;
    char    dname[NC_MAX_NAME + 1];

    if (strcmp(c->lai_fname, "*NOT SET*") == 0) {
        return (read_nc_var(met_ncid, "LAI", n, TRUE));
    }

    NC_CHECK(nc_open(c->lai_fname, NC_NOWRITE, &ncid), c->lai_fname);
    NC_CHECK(nc_inq_varid(ncid, c->lai_var, &varid), c->lai_var);
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims), c->lai_var);
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids), c->lai_var);
    for (d = 0; d < ndims; d++) {
        NC_CHECK(nc_inq_dim(ncid, dimids[d], dname, &len), c->lai_var);
        start[d] = 0;
        count[d] = 1;
        if (strcmp(dname, "time") == 0) {
            idim_time = d;
            if (len != n) {
                // not on the met timesteps: a daily file
                lai_day = read_lai_days(ncid, len);
                nlai = len;
            }
            count[d] = nlai;
        } else if (strcmp(dname, "pft") == 0) {
            start[d] = (size_t)c->lai_pft_index;
        }
    }
    if (idim_time < 0) {
        fprintf(stderr, "LAI variable %s has no time dimension\n", c->lai_var);
        exit(EXIT_FAILURE);
    }
    daily = alloc_array(nlai, "LAI");
    NC_CHECK(nc_get_vara_double(ncid, varid, start, count, daily),
             c->lai_var);
    nc_close(ncid);
    if (lai_day == NULL) {
        out = daily;
    } else {
        out = alloc_array(n, "LAI");
        for (i = 0; i < n; i++) {
            k = met_day[i] - lai_day[0];
            if (k < 0 || k >= (long)nlai) {
                fprintf(stderr, "Daily LAI file doesn't cover met timestep "
                        "%zu\n", i);
                exit(EXIT_FAILURE);
            }
            out[i] = daily[k];
        }
        free(daily);
        free(lai_day);
    }
    for (i = 0; i < n; i++) {
        if (out[i] < 0.0 || isnan(out[i])) {
            fprintf(stderr, "Bad LAI value at timestep %zu\n", i);
            exit(EXIT_FAILURE);
        }
    }

    return (out);
}

static double qair_to_vpd(double qair, double tair, double press) {
    /* VPD (kPa) from specific humidity, tair (deg C), press (Pa) */
    double es = 100.0 * 6.112 * exp((17.67 * tair) / (243.5 + tair));
    double ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair);
    return (MAX(0.05, (es - ea) * PA_2_KPA));
}

static double *alloc_array(size_t n, const char *name) {
    double *a = calloc(n, sizeof(double));
    if (a == NULL) {
        fprintf(stderr, "Error allocating space for %s array\n", name);
        exit(EXIT_FAILURE);
    }
    return (a);
}

void read_met_data_netcdf(char **argv, control *c, met_arrays *ma,
                          params *p) {

    int     ncid, dimid, varid, y0, mo0, d0, h0, mi0, s0;
    size_t  ntime, i, i0, i1, n, j, k, nday, step, spd;
    double *time, *tair, *sw, *precip, *qair, *psurf, *wind, *co2, *lai, *lwdown;
    double  dt, t0_days, current_yr, per_step, par;
    char    units[NC_MAX_NAME + 1];
    int    *yr_of, *doy_of;
    long   *met_day;

    NC_CHECK(nc_open(c->met_fname, NC_NOWRITE, &ncid), c->met_fname);
    NC_CHECK(nc_inq_dimid(ncid, "time", &dimid), "time dimension");
    NC_CHECK(nc_inq_dimlen(ncid, dimid, &ntime), "time dimension");
    NC_CHECK(nc_inq_varid(ncid, "time", &varid), "time variable");

    memset(units, 0, sizeof(units));
    NC_CHECK(nc_get_att_text(ncid, varid, "units", units), "time units");
    if (sscanf(units, "seconds since %d-%d-%d %d:%d:%d", &y0, &mo0, &d0, &h0,
               &mi0, &s0) != 6) {
        fprintf(stderr, "Can't parse time units '%s' (expected 'seconds "
                "since YYYY-MM-DD hh:mm:ss')\n", units);
        exit(EXIT_FAILURE);
    }
    time = alloc_array(ntime, "time");
    NC_CHECK(nc_get_var_double(ncid, varid, time), "time");
    dt = time[1] - time[0];
    spd = (size_t)(86400.0 / dt + 0.5);
    if (c->sub_daily && spd != (size_t)c->num_hlf_hrs) {
        fprintf(stderr, "The sub-daily model needs 30 min forcing, the "
                "netCDF timestep is %.0f s\n", dt);
        exit(EXIT_FAILURE);
    }

    /* date of each timestep (timestamps mark the start of the interval) */
    t0_days = (double)days_from_civil(y0, mo0, d0) +
              (h0 * 3600.0 + mi0 * 60.0 + s0) / 86400.0;
    yr_of = malloc(ntime * sizeof(int));
    doy_of = malloc(ntime * sizeof(int));
    met_day = malloc(ntime * sizeof(long));
    for (i = 0; i < ntime; i++) {
        met_day[i] = (long)floor(t0_days + time[i] / 86400.0 + 1E-9);
        civil_from_days(met_day[i], &yr_of[i], &doy_of[i]);
    }

    /* optional subset of whole years */
    i0 = 0;
    i1 = ntime;
    if (c->met_start_year > 0) {
        while (i0 < ntime && yr_of[i0] < c->met_start_year) i0++;
    }
    if (c->met_end_year > 0) {
        while (i1 > i0 && yr_of[i1 - 1] > c->met_end_year) i1--;
    }
    n = i1 - i0;
    if (n == 0 || n % spd != 0) {
        fprintf(stderr, "netCDF forcing must contain whole days (%zu steps)\n",
                n);
        exit(EXIT_FAILURE);
    }

    tair = read_nc_var(ncid, "Tair", ntime, TRUE);
    sw = read_nc_var(ncid, "SWdown", ntime, TRUE);
    precip = read_nc_var(ncid, "Precip", ntime, TRUE);
    qair = read_nc_var(ncid, "Qair", ntime, TRUE);
    psurf = read_nc_var(ncid, "Psurf", ntime, TRUE);
    wind = read_nc_var(ncid, "Wind", ntime, TRUE);
    co2 = read_nc_var(ncid, "CO2air", ntime, FALSE);
    lwdown = c->sub_daily ? read_nc_var(ncid, "LWdown", ntime, FALSE) : NULL;
    lai = c->prescribed_lai ? read_nc_lai(c, ntime, ncid, met_day) : NULL;
    free(met_day);
    nc_close(ncid);

    if (co2 == NULL && p->nc_co2 < 0.0) {
        fprintf(stderr, "netCDF met file has no CO2air, set nc_co2\n");
        exit(EXIT_FAILURE);
    }

    nday = n / spd;
    c->num_years = 0;
    current_yr = -999.9;
    ma->lai = NULL;
    ma->lwdown = NULL;

    if (c->sub_daily) {
        c->total_num_days = nday;
        ma->year = alloc_array(n, "year");
        ma->doy = alloc_array(n, "doy");
        ma->rain = alloc_array(n, "rain");
        ma->par = alloc_array(n, "par");
        ma->tair = alloc_array(n, "tair");
        ma->tsoil = alloc_array(n, "tsoil");
        ma->vpd = alloc_array(n, "vpd");
        ma->co2 = alloc_array(n, "co2");
        ma->ndep = alloc_array(n, "ndep");
        ma->nfix = alloc_array(n, "nfix");
        ma->wind = alloc_array(n, "wind");
        ma->press = alloc_array(n, "press");
        if (lai != NULL) {
            ma->lai = alloc_array(n, "lai");
        }
        if (lwdown != NULL) {
            ma->lwdown = alloc_array(n, "lwdown");
        }

        per_step = dt / (NDAYS_IN_YR * 86400.0);
        for (k = 0; k < nday; k++) {
            double tmean = 0.0;
            for (step = 0; step < spd; step++) {
                tmean += tair[i0 + k * spd + step] - DEG_TO_KELVIN;
            }
            tmean /= (double)spd;

            for (step = 0; step < spd; step++) {
                i = i0 + k * spd + step;
                j = k * spd + step;
                ma->year[j] = yr_of[i];
                ma->doy[j] = doy_of[i];
                ma->rain[j] = MAX(0.0, precip[i] * dt);         /* mm */
                ma->par[j] = MAX(0.0, sw[i] * SW_2_PAR);  /* umol m-2 s-1 */
                ma->tair[j] = tair[i] - DEG_TO_KELVIN;
                ma->tsoil[j] = tmean;
                ma->vpd[j] = qair_to_vpd(qair[i], ma->tair[j], psurf[i]);
                ma->co2[j] = co2 != NULL ? co2[i] : p->nc_co2;
                ma->ndep[j] = p->nc_ndep * per_step;
                ma->nfix[j] = p->nc_nfix * per_step;
                ma->wind[j] = wind[i];
                ma->press[j] = psurf[i] * PA_2_KPA;
                if (lai != NULL) {
                    ma->lai[j] = lai[i];
                }
                if (lwdown != NULL) {
                    ma->lwdown[j] = lwdown[i];
                }
                if (current_yr != ma->year[j]) {
                    c->num_years++;
                    current_yr = ma->year[j];
                }
            }
        }
    } else {
        c->total_num_days = nday;
        ma->year = alloc_array(nday, "year");
        ma->prjday = alloc_array(nday, "prjday");
        ma->tair = alloc_array(nday, "tair");
        ma->rain = alloc_array(nday, "rain");
        ma->tsoil = alloc_array(nday, "tsoil");
        ma->tam = alloc_array(nday, "tam");
        ma->tpm = alloc_array(nday, "tpm");
        ma->tmin = alloc_array(nday, "tmin");
        ma->tmax = alloc_array(nday, "tmax");
        ma->tday = alloc_array(nday, "tday");
        ma->vpd_am = alloc_array(nday, "vpd_am");
        ma->vpd_pm = alloc_array(nday, "vpd_pm");
        ma->co2 = alloc_array(nday, "co2");
        ma->ndep = alloc_array(nday, "ndep");
        ma->nfix = alloc_array(nday, "nfix");
        ma->wind = alloc_array(nday, "wind");
        ma->press = alloc_array(nday, "press");
        ma->wind_am = alloc_array(nday, "wind_am");
        ma->wind_pm = alloc_array(nday, "wind_pm");
        ma->par = alloc_array(nday, "par");
        ma->par_am = alloc_array(nday, "par_am");
        ma->par_pm = alloc_array(nday, "par_pm");
        if (lai != NULL) {
            ma->lai = alloc_array(nday, "lai");
        }

        for (k = 0; k < nday; k++) {
            /* sums: daylight, am, pm; T, wind, press, qair, co2; counts */
            double s_t = 0, s_w = 0, s_p = 0, s_c = 0, n_l = 0;
            double a_t = 0, a_w = 0, a_p = 0, a_q = 0, a_par = 0, n_a = 0;
            double p_t = 0, p_w = 0, p_p = 0, p_q = 0, p_par = 0, n_p = 0;
            double t24 = 0, tmin = 999.9, tmax = -999.9, rain = 0, s_lai = 0;
            double tc, conv = UMOL_2_JOL * J_TO_MJ * dt; /* umol s-1 -> MJ */

            for (step = 0; step < spd; step++) {
                i = i0 + k * spd + step;
                tc = tair[i] - DEG_TO_KELVIN;
                par = MAX(0.0, sw[i] * SW_2_PAR);
                t24 += tc;
                tmin = MIN(tmin, tc);
                tmax = MAX(tmax, tc);
                rain += MAX(0.0, precip[i] * dt);
                if (lai != NULL) {
                    s_lai += lai[i];
                }
                if (par >= 5.0) {
                    s_t += tc;
                    s_w += wind[i];
                    s_p += psurf[i];
                    s_c += co2 != NULL ? co2[i] : p->nc_co2;
                    n_l += 1.0;
                    if (step < spd / 2) {
                        a_t += tc; a_w += wind[i]; a_p += psurf[i];
                        a_q += qair[i]; a_par += par * conv; n_a += 1.0;
                    } else {
                        p_t += tc; p_w += wind[i]; p_p += psurf[i];
                        p_q += qair[i]; p_par += par * conv; n_p += 1.0;
                    }
                }
            }
            i = i0 + k * spd;
            ma->year[k] = yr_of[i];
            ma->prjday[k] = doy_of[i];
            ma->tday[k] = t24 / (double)spd;
            ma->tsoil[k] = ma->tday[k];
            ma->tmin[k] = tmin;
            ma->tmax[k] = tmax;
            ma->rain[k] = rain;
            ma->ndep[k] = p->nc_ndep / NDAYS_IN_YR;
            ma->nfix[k] = p->nc_nfix / NDAYS_IN_YR;
            if (lai != NULL) {
                ma->lai[k] = s_lai / (double)spd;
            }

            /* polar night etc: fall back on the 24 h values */
            ma->tair[k] = n_l > 0 ? s_t / n_l : ma->tday[k];
            ma->wind[k] = MAX(0.1, n_l > 0 ? s_w / n_l : wind[i]);
            ma->press[k] = (n_l > 0 ? s_p / n_l : psurf[i]) * PA_2_KPA;
            ma->co2[k] = n_l > 0 ? s_c / n_l : (co2 != NULL ? co2[i]
                                                            : p->nc_co2);
            ma->tam[k] = n_a > 0 ? a_t / n_a : ma->tair[k];
            ma->tpm[k] = n_p > 0 ? p_t / n_p : ma->tair[k];
            ma->wind_am[k] = MAX(0.1, n_a > 0 ? a_w / n_a : ma->wind[k]);
            ma->wind_pm[k] = MAX(0.1, n_p > 0 ? p_w / n_p : ma->wind[k]);
            ma->vpd_am[k] = n_a > 0 ? qair_to_vpd(a_q / n_a, ma->tam[k],
                                                  a_p / n_a) : 0.05;
            ma->vpd_pm[k] = n_p > 0 ? qair_to_vpd(p_q / n_p, ma->tpm[k],
                                                  p_p / n_p) : 0.05;
            ma->par_am[k] = a_par;
            ma->par_pm[k] = p_par;

            if (current_yr != ma->year[k]) {
                c->num_years++;
                current_yr = ma->year[k];
            }
        }
    }

    free(time); free(tair); free(sw); free(precip); free(qair); free(psurf);
    free(wind); free(co2); free(lai); free(lwdown); free(yr_of); free(doy_of);
    (void)argv;

    return;
}

#else

void read_met_data_netcdf(char **argv, control *c, met_arrays *ma,
                          params *p) {
    fprintf(stderr, "%s: this GDAY was built without netCDF support, "
            "rebuild with netCDF (see the Makefile) to read %s\n",
            argv[0], c->met_fname);
    exit(EXIT_FAILURE);
}

#endif
