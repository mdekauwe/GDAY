#include "read_met_file.h"

void read_daily_met_data(char **argv, control *c, met_arrays *ma)
{
    FILE  *fp;
    char   line[STRING_LENGTH];
    int    file_len = 0;
    int    i = 0;
    int    nvars = 21;
    int    skipped_lines = 0;
    double current_yr = -999.9;

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open daily Met file %s for read\n",
                c->met_fname);
		exit(EXIT_FAILURE);
	 }

    /* work out how big the file is */
    file_len = 0;
    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#')
            continue;
        file_len++;
    }
    rewind(fp);
    c->total_num_days = file_len;

    /* allocate memory for meteorological arrays */
    if ((ma->year = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for year array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->prjday = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for prjday array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tair = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tair array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->rain = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for rain array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tsoil = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tsoil array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tam = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tam array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tpm = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tpm array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tmin = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tmin array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tmax = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tmax array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tday = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tday array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->vpd_am = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->vpd_pm = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd_pm array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->co2 = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for co2 array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->ndep = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for ndep array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->nfix = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for nfix array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->wind = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->press = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for press array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->wind_am = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->wind_pm = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind_pm array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par_am = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par_pm = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par_pm array\n");
		exit(EXIT_FAILURE);
    }


    i = 0;
    c->num_years = 0;
    skipped_lines = 0;
    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#') {
            skipped_lines++;
            continue;
        }

        if (sscanf(line, "%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,\
                          %lf,%lf",\
                          &(ma->year[i]), &(ma->prjday[i]), \
                          &(ma->tair[i]), &(ma->rain[i]), &(ma->tsoil[i]), \
                          &(ma->tam[i]), &(ma->tpm[i]), &(ma->tmin[i]), \
                          &(ma->tmax[i]), &(ma->tday[i]), &(ma->vpd_am[i]), \
                          &(ma->vpd_pm[i]), &(ma->co2[i]), &(ma->ndep[i]), \
                          &(ma->nfix[i]),  &(ma->wind[i]), &(ma->press[i]), \
                          &(ma->wind_am[i]), &(ma->wind_pm[i]), \
                          &(ma->par_am[i]), &(ma->par_pm[i])) != nvars) {
            fprintf(stderr, "%s: badly formatted input in met file on line %d %d\n", \
                    *argv, (int)i+1+skipped_lines, nvars);
            exit(EXIT_FAILURE);
        }

        /* Build an array of the unique years as we loop over the input file */
        if (current_yr != ma->year[i]) {
            c->num_years++;
            current_yr = ma->year[i];
        }
        i++;
    }
    fclose(fp);
    return;
}

void read_subdaily_met_data(char **argv, control *c, met_arrays *ma)
{
    FILE  *fp;
    char   line[STRING_LENGTH];
    int    i = 0;
    int    nvars = 13;
    int    skipped_lines = 0;
    double current_yr, temp_HOD;
    long   file_len;

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open sub-daily Met file %s for read\n",
                c->met_fname);
		exit(EXIT_FAILURE);
	 }

    /* work out how big the file is */
    file_len = 0;
    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#')
            continue;
        file_len++;
    }
    rewind(fp);

    /* output is daily, so correct for n_timesteps */
    c->total_num_days = file_len / 48;

    /* allocate memory for meteorological arrays */
    if ((ma->year = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for year array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->doy = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for doy array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->rain = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for rain array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tair = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tair array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tsoil = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tsoil array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->vpd = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->co2 = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for co2 array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->ndep = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for ndep array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->nfix = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for nfix array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->wind = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->press = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for press array\n");
		exit(EXIT_FAILURE);
    }

    current_yr = ma->year[0];

    i = 0;
    c->num_years = 0;
    skipped_lines = 0;

    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#') {
            skipped_lines++;
            continue;
        }

        if (sscanf(line, "%lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf", \
                          &(ma->year[i]), &(ma->doy[i]), &temp_HOD, \
                          &(ma->rain[i]), &(ma->par[i]), &(ma->tair[i]), \
                          &(ma->tsoil[i]), &(ma->vpd[i]), &(ma->co2[i]), \
                          &(ma->ndep[i]), &(ma->nfix[i]), &(ma->wind[i]), \
                          &(ma->press[i])) != nvars) {
            fprintf(stderr, "%s: badly formatted input in met file on line %d %d\n", \
                    *argv, (int)i+1+skipped_lines, nvars);
            exit(EXIT_FAILURE);
        }

        /* Build an array of the unique years as we loop over the input file */
        if (current_yr != ma->year[i]) {
            c->num_years++;
            current_yr = ma->year[i];
        }
        i++;
    }

    fclose(fp);
    return;
}


void read_daily_met_data_binary(char **argv, control *c, met_arrays *ma)
{
    return;
}

void read_subdaily_met_data_binary(char **argv, control *c, met_arrays *ma,
                                   params *p, state *s)
{
    FILE  *fp;
    float *data = NULL;
    int    i = 0, hod, doy_idx;
    double current_yr;
    long   cnt;
    double tmin, tmax, vph09, vph15, sw, rain, co2, tsoil;
    double vph09_tomorrow, vph15_yesterday;

    // Need to unpack everything into a 48 hour array
    float parx[NHRS], vphx[NHRS], tairx[NHRS], rainx[NHRS];

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open Met file %s for read\n",
                c->met_fname);
		exit(EXIT_FAILURE);
    }

    if ((data = (float *)calloc(c->nrows*c->ncols, sizeof(float))) == NULL) {
		fprintf(stderr,"%s: error allocating met data for read\n",argv[0]);
		exit(EXIT_FAILURE);
	}

    if (fread(data, sizeof(float), c->nrows*c->ncols, fp) != c->nrows*c->ncols){
		fprintf(stderr,"Error in reading met binary file: %s\n", c->met_fname);
		exit(EXIT_FAILURE);
	}

    c->total_num_days = c->nrows;

    /* allocate memory for meteorological arrays */
    if ((ma->year = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for year array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->doy = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for doy array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->rain = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for rain array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->par = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tair = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tair array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->tsoil = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tsoil array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->vpd = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->co2 = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for co2 array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->ndep = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for ndep array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->nfix = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for nfix array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->wind = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->press = (double *)calloc(c->nrows * 48, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for press array\n");
		exit(EXIT_FAILURE);
    }

    /* Unpack the met chunk */
    cnt = 0;
    c->num_years = 1;
    current_yr = data[0];
    for (i = 0; i < c->nrows * c->ncols; i += c->ncols) {

        ma->year[cnt] = data[i];
        ma->doy[cnt] = data[i+1];
        doy_idx = (int)ma->doy[cnt] - 1;
        tmin = data[i+2];
        tmax = data[i+3];
        vph09 = data[i+4];
        vph15 = data[i+5];
        vph09_tomorrow = data[i+6];
        vph15_yesterday = data[i+7];
        sw = data[i+8];
        rain = data[i+9];
        co2 = data[i+10];

        estimate_dirunal_par(p->latitude, p->longitude, ma->doy[cnt], sw,
                             &(parx[0]));
        estimate_diurnal_vph(vph09, vph15, vph09_tomorrow, vph15_yesterday,
                             &(vphx[0]));
        disaggregate_rainfall(rain, &(rainx[0]));
        estimate_diurnal_temp(tmin, tmax, s->day_length[doy_idx],
                              &(tairx[0]));

        for (hod = 0; hod < NHRS; hod++) {
            tsoil += tairx[hod];
        }
        tsoil /= (float)NHRS;

        // Unpacking met data ...
        for (hod = 0; hod < NHRS; hod++) {

            ma->rain[i] = rainx[hod];
            ma->par[i] = parx[hod];
            ma->tair[i] = tairx[hod];
            ma->tsoil[i] = tsoil;
            ma->vpd[i] = calc_vpd(tairx[hod], vphx[hod]);
            ma->co2[i] = co2;
            ma->ndep[i] = -9999.9;
            ma->nfix[i] = -9999.9;
            ma->wind[i] = 3.0;    /* Haverd et al. 2012 */
            ma->press[i] = 100.0;
        }

        /* Build an array of the unique years as we loop over the input file */
        if (current_yr != ma->year[cnt]) {
            c->num_years++;
            current_yr = ma->year[cnt];
        }
        cnt++;

    }
    exit(1);


    /* tidy up */
    free(data);
    fclose(fp);

    return;
}

void estimate_dirunal_par(float lat, float lon, int doy, float sw_rad_day,
                          float *par) {
    /*
        Calculate daily course of incident PAR from daily totals using routine
        from MAESTRA
    */
    int   i;
    float cos_zenith[NHRS];
    float tau = 0.76;            /* Transmissivity of atmosphere */
    float direct_frac, diffuse_frac;
    float cos_bm[NHRS], cos_df[NHRS], sum_bm, sum_df;
    float zenith, rddf, rdbm, par_day, beam_rad, diffuse_rad;

    /* MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> umol m-2 s-1 -> MJ m-2 d-1 */
    par_day = sw_rad_day * MJ_TO_J * DAY_2_SECX * SW_2_PAR * \
              UMOL_TO_J * J_TO_MJ * SEC_2_DAYX;

    calculate_solar_geometryx(doy, lat, lon, &(cos_zenith[0]));
    diffuse_frac = spittersx(doy, par_day, cos_zenith);
    direct_frac = 1.0 - diffuse_frac;

    /* daily total beam PAR (MJ m-2 d-1) */
    beam_rad = par_day * direct_frac;

    /* daily total diffuse PAR (MJ m-2 d-1) */
    diffuse_rad = par_day * diffuse_frac;

    sum_bm = 0.0;
    sum_df = 0.0;
    for (i = 0; i < NHRS; i++) {
        cos_bm[i] = 0.0;
        cos_df[i] = 0.0;

        if (cos_zenith[i] > 0.0) {
            zenith = acos(cos_zenith[i]);

            /* set FBM = 0.0 for ZEN > 80 degrees */
            if (zenith < (80.0 * M_PI / 180.0)) {
                cos_bm[i] = cos_zenith[i] * pow(tau, (1.0 / cos_zenith[i]));
            } else {
                cos_bm[i] = 0.0;
            }
            cos_df[i] = cos_zenith[i];
            sum_bm += cos_bm[i];
            sum_df += cos_df[i];
        }
    }

    for (i = 0; i < NHRS; i++) {

        if (sum_bm > 0.0) {
            rdbm = beam_rad * cos_bm[i] / sum_bm;
        } else {
            rdbm = 0.0;
        }

        if (sum_df > 0.0) {
            rddf = diffuse_rad * cos_df[i] / sum_df;
        } else {
            rddf = 0.0;
        }

        /* MJ m-2 30min-1 -> J m-2 s-1 -> umol m-2 s-1 */
        *(par+i) = (rddf + rdbm) * MJ_TO_J * J_TO_UMOL * HLFHR_2_SEC;
    }

    return;
}

void estimate_diurnal_vph(float vph09, float vph15, float vph09_next,
                          float vph15_prev, float *vph) {
    /*
    Interpolate VPH between 9am and 3pm values to generate diurnal VPD
    following the method of Haverd et al. This seems reasonable, vapour pressure
    plotted aginst time of day often does not reveal consistent patterns, with
    small fluctuations (see Kimball and Bellamy, 1986).
    Reference:
    ---------
    * Haverd et al. (2013) Multiple observation types reduce uncertainty in
      Australia's terrestrial carbon and water cycles. Biogeosciences, 10,
      2011-2040.
    */
    /* number of hours gap, i.e. 3pm to 9am the next day */
    float gap = 18.0;
    float hour;
    int    i;

    for (i = 1; i < NHRS+1; i++) {
        /* first zero values */
        *(vph+i) = 0.0;

        hour = (float)i / 2.0;

        if (hour <= 9.0) {
           *(vph+(i-1)) = vph15_prev + (vph09 - vph15_prev) * (9.0 + hour) / gap;
       } else if (hour > 9.0 && hour <= 15.0) {
           *(vph+(i-1)) = vph09 + (vph15 - vph09) * (hour - 9.0) / (15.0 - 9.0);
        } else if (hour > 15.0) {
            *(vph+(i-1)) =  vph15 + (vph09_next - vph15) * (hour - 15.0) / gap;
        }
    }

    return;
}

void disaggregate_rainfall(float rain_day, float *rain) {
    /*
    Assign daily PPT total to hours of the day, following MAESTRA, which follows
    algorithm from GRAECO (model of D. Loustau).
    Reference:
    * Loustau, D., F. Pluviaud, A. Bosc, A. Porte, P. Berbigier, M. Deque
      and V. Perarnaud. 2001. Impact of a regional 2 x CO2 climate scenario
      on the water balance, carbon balance and primary production
      of maritime pine in southwestern France. In Models for the Sustainable
      Management of Plantation Forests. Ed. M. Tome. European
      Cultivated Forest Inst., EFI Proc. No. 41D, Bordeaux, pp 45-58.
    */
    int   i, j, hour_index, num_hrs_with_rain;
    float rate;

    /* zero everything before we start */
    for (i = 0; i < NHRS; i++) {
        *(rain+i) = 0.0;
    }

    if (rain_day <= 2.0) {
        /* All rain falls in one hour for light storms (<2 mm) */
        hour_index = rand_int(0, 47);
        *(rain+hour_index) = rain_day;

    } else if (rain_day > 46.0) {
        /* All rain falls in 24 hours for storms >46 mm */
        for (i = 0; i < NHRS; i++) {
            *(rain+i) = rain_day / (float)NHRS;
        }

    } else {
        /*
        ** Aim if for all rain to fall at 2mm/hour at a random time of the day.
        ** If we generate the same random number, then we increase rainfall
        ** for this hour
        */
        num_hrs_with_rain = (int)(rain_day / 2.0);
        rate = rain_day / (float)num_hrs_with_rain;

        for (j = 0; j < num_hrs_with_rain; j++) {
            hour_index = rand_int(0, 47);
            *(rain+hour_index) += rate;
        }
    }

    return;
}

void estimate_diurnal_temp(float tmin, float tmax, float day_length,
                           float *tair) {
    /*
    Calculate diurnal temperature following Parton and Logan
    the day is divided into two segments and using a truncated sine wave
    in the daylight and an exponential decrease in temperature
    at night.
    TO DO:
    - Hours between 00:00 and sunrise should be modelled using the previous
      days information.
    References:
    ----------
    * Parton and Logan (1981) A model for dirunal variation in soil and
       air temperature. Agricultural Meteorology, 23, 205--216.
    * Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197.
    */
    /* 1.5 m air temperature values from Parton and Logan, table 1 */
    float a = 1.86;
    float b = 2.2;     /* nighttime coeffcient */
    float c = -0.17;   /* lag of the min temp from the time of runrise */

    float night_length = 24.0 - day_length;
    float sunrise = 12.0 - day_length / 2.0 + c;
    float sunset = 12.0 + day_length / 2.0;
    float m, n, d, tset, hour;
    int   i;

    /* temperature at sunset */
    m = sunset - sunrise + c;
    tset = (tmax - tmin) * sin(M_PI * m / (day_length + 2.0 * a)) + tmin;


    for (i = 1; i < NHRS+1; i++) {

        hour = (float)i / 2.0;

        /* hour - time of the minimum temperature (accounting for lag time) */
        m = hour - sunrise + c;
        if (hour >= sunrise && hour <= sunset) {
            *(tair+(i-1)) = tmin + (tmax - tmin) * \
                        sin((M_PI * m) / (day_length + 2.0 * a));
        } else {
            if (hour > sunset) {
                n = hour - sunset;
            } else if (hour < sunrise) {
                n = (24.0 + hour) - sunset;
            }

            d = (tset - tmin) / (exp(b) - 1.0);

            /* includes missing displacement to allow T to reach Tmin, this
            ** removes a discontinuity in the original Parton and Logan eqn.
            ** See Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197
            **/
            *(tair+(i-1)) = (tmin -d) + (tset - tmin - d) * \
                        exp(-b * n / (night_length + c));
        }
    }

    return;
}

float calc_vpd(float temp, float ea) {

    /*
    Empirical equation following Tetens (1930), using the form of Murray
    (1967) as given in Montieth and Unsworth (1990), pg. 10.
    Parameters:
    -----------
    temp : float
        air temperature (deg C)
    References:
    ----------
    * Monteith JL & Unsworth MH (1990) Principles of environmental physics.
    */

    float Tk, A, T_star, T_dash, es_T_star, esat, vpd;

    Tk = temp + DEG_TO_KELVIN;

    A = 17.27;
    T_star = 273.0; /* K */
    T_dash = 36.0;  /* K */
    es_T_star = 0.611; /* kPa */

    /* saturation vapor pressure (kPa)
       Values of saturation vapour pressure from the Tetens formula are
       within 1 Pa of the exact values.
    */
    esat = es_T_star * exp(A * (Tk - T_star) / (Tk - T_dash));

    /* VPD is the saturated vapour pressure - the actual vapour pressure */
    vpd = MAX(0.05, esat - (ea * hPa_2_kPa));

    return vpd;
}

int rand_int(unsigned int min, unsigned int max) {

    int   value;
    float scaled;

    scaled = (float)rand() / RAND_MAX;
    value = (int)(max - min + 1) * scaled + min;

    return (value);
}

void calculate_solar_geometryx(int doy, float latitude, float longitude,
                              float *cos_zenith) {
    /*
    The solar zenith angle is the angle between the zenith and the centre
    of the sun's disc. The solar elevation angle is the altitude of the
    sun, the angle between the horizon and the centre of the sun's disc.
    Since these two angles are complementary, the cosine of either one of
    them equals the sine of the other, i.e. cos theta = sin beta. I will
    use cos_zen throughout code for simplicity.
    Arguments:
    ----------
    doy : float
        day of year
    latitude : float
        latitude (degrees)
    longitude : float
        longitude (degrees)
    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    */
    int   i;
    float rdec, et, t0, h, gamma, rlat, sin_beta;
    float hod;

    for (i = 1; i < NHRS+1; i++) {

        /* need to convert 30 min data, 0-47 to 0-23.5 */
        hod = i / 2.0;

        gamma = day_anglex(doy);
        rdec = calculate_solar_declinationx(doy, gamma);
        et = calculate_eqn_of_timex(gamma);
        t0 = calculate_solar_noonx(et, longitude);
        h = calculate_hour_anglex(hod, t0);
        rlat = latitude * M_PI / 180.0;

        /* A13 - De Pury & Farquhar */
        sin_beta = sin(rlat) * sin(rdec) + cos(rlat) * cos(rdec) * cos(h);
        /* The same thing, going to use throughout */
        *(cos_zenith+(i-1)) = sin_beta;
        if (*(cos_zenith+(i-1)) > 1.0) {
            *(cos_zenith+(i-1)) = 1.0;
        } else if (cos_zenith[i-1] < 0.0) {
            *(cos_zenith+(i-1)) = 0.0;
        }
        /*zenith = 180.0 / M_PI * acos(cos_zenith[i-1]);
        elevation = 90.0 - zenith;*/
    }
    return;
}

float day_anglex(int doy) {
    /* Calculation of day angle - De Pury & Farquhar, '97: eqn A18
    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    Returns:
    ---------
    gamma - day angle in radians.
    */

    return (2.0 * M_PI * ((float)doy - 1.0) / 365.0);
}

float calculate_solar_declinationx(int doy, float gamma) {
    /*
    Solar Declination Angle is a function of day of year and is indepenent
    of location, varying between 23deg45' to -23deg45'
    Arguments:
    ----------
    doy : int
        day of year, 1=jan 1
    gamma : float
        fractional year (radians)
    Returns:
    --------
    dec: float
        Solar Declination Angle [radians]
    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    */
    float decl;

    /* declination (radians) */
    /*decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
           0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
           0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/


    /* (radians) A14 - De Pury & Farquhar  */
    decl = -23.4 * (M_PI / 180.) * cos(2.0 * M_PI * ((float)doy + 10.) / 365.);

    return (decl);

}

float calculate_eqn_of_timex(float gamma) {
    /* Equation of time - correction for the difference btw solar time
    and the clock time.
    Arguments:
    ----------
    doy : int
        day of year
    gamma : float
        fractional year (radians)
    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
      biophysics. Pg 169.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
      "The Equation of Time", Monthly Notices of the Royal Astronomical
      Society 238: 1529–1535
    */
    float et;

    /* radians */
    et = 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -\
         0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma);

    /* radians to minutes */
    et *= 229.18;

    /* radians to hours */
    /*et *= 24.0 / (2.0 * M_PI);*/

    /* minutes - de Pury and Farquhar, 1997 - A17 */
    /*et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
          cos(2.0 * gamma) - 9.731  * sin(gamma));*/

    return (et);
}

float calculate_solar_noonx(float et, float longitude) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A16
    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    Returns:
    ---------
    t0 - solar noon (hours).
    */
    float t0, Ls;

    /* all international standard meridians are multiples of 15deg east/west of
       greenwich */
    Ls = round_to_valuex(longitude, 15.);
    t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0;

    return (t0);
}

float calculate_hour_anglex(float t, float t0) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A15
    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    Returns:
    ---------
    h - hour angle (radians).
    */
    return (M_PI * (t - t0) / 12.0);

}


float spittersx(int doy, float par, float *cos_zenith) {
    /*
    Spitters algorithm to estimate the diffuse component from the total daily
    incident radiation.
    NB. Eqns. 2a-d, not 20a-d
    Parameters:
    ----------
    doy : int
        day of year
    par : float
        daily total photosynthetically active radiation (MJ m-2 d-1)
    cos_zenith : float
        cosine of zenith angle (radians)
    Returns:
    -------
    diffuse : float
        diffuse component of incoming radiation
    References:
    ----------
    * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
      Separating the diffuse and direct component of global radiation and its
      implications for modeling canopy photosynthesis. Part I. Components of
      incoming radiation. Agricultural Forest Meteorol., 38:217-229.
    */

    /* Fraction of global radiation that is PAR */
    float fpar = 0.5;
    float conv = SEC_2_HFHR * J_TO_MJ;
    float S0, tau, diffuse_frac;
    int   i;


    /* Calculate extra-terrestrial radiation */
    S0 = 0.0;
    for (i = 1; i < NHRS+1; i++) {
        S0 += calc_extra_terrestrial_radx(doy, *(cos_zenith+(i-1))) * conv;
    }

    /* atmospheric transmisivity */
    tau = (par / fpar) / S0;

    /* Spitter's formula (Eqns. 2a-d) */
    if (tau < 0.07) {
        diffuse_frac = 1.0;
    } else if (tau < 0.35) {
        diffuse_frac = 1.0 - 2.3 * (tau - 0.07) * (tau - 0.07);
    } else if (tau < 0.75) {
        diffuse_frac = 1.33 - 1.46 * tau;
    } else {
        diffuse_frac = 0.23;
    }

    return (diffuse_frac);
}

float calc_extra_terrestrial_radx(int doy, float cos_zenith) {
    /* Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the earths
    orbit.
    Using formula from Spitters not Leuning!
    Arguments:
    ----------
    doy : double
        day of year
    cos_zenith : double
        cosine of zenith angle (radians)
    Returns:
    --------
    So : float
        solar radiation normal to the sun's bean outside the Earth's atmosphere
        (J m-2 s-1)
    Reference:
    ----------
    * Spitters et al. (1986) AFM, 38, 217-229, equation 1.
    */

    float So, Sc;

    /* Solar constant (J m-2 s-1) */
    Sc = 1370.0;

    if (cos_zenith > 0.0) {
        /*
        ** remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
        ** sin(x) = cos(90-x) and cos(x) = sin(90-x).
        */
        So = Sc * (1.0 + 0.033 * cos((float)doy / 365.0 * 2.0 * M_PI)) *\
                cos_zenith;
    } else {
        So = 0.0;
    }

    return (So);

}

float round_to_valuex(float number, float roundto) {
    return (round(number / roundto) * roundto);
}
