#include "read_met_file.h"

void read_daily_met_data(char **argv, control *c, met_arrays *ma)
{
    FILE  *fp;
    char   line[STRING_LENGTH];
    int    file_len = 0;
    int    i = 0;
    int    nvars = 20;
    int    skipped_lines = 0;
    double current_yr;

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open Met file %s for read\n",
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
                          %lf,%lf,%lf", \
                          &(ma->year[i]), &(ma->prjday[i]), \
                          &(ma->tair[i]), &(ma->rain[i]), &(ma->tsoil[i]), \
                          &(ma->tam[i]), &(ma->tpm[i]), &(ma->tmin[i]), \
                          &(ma->tmax[i]), &(ma->tday[i]), &(ma->vpd_am[i]), \
                          &(ma->vpd_pm[i]), &(ma->co2[i]), &(ma->ndep[i]), \
                          &(ma->wind[i]), &(ma->press[i]), &(ma->wind_am[i]), \
                          &(ma->wind_pm[i]), &(ma->par_am[i]), \
                          &(ma->par_pm[i])) != nvars) {
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
    int    nvars = 12;
    int    skipped_lines = 0;
    double current_yr, temp_HOD;
    long   file_len;
    double hod_dummy;

    if ((fp = fopen(c->met_fname, "r")) == NULL) {
		fprintf(stderr, "Error: couldn't open Met file %s for read\n",
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
                          %lf,%lf,%lf", \
                          &(ma->year[i]), &(ma->doy[i]), &temp_HOD, \
                          &(ma->rain[i]), &(ma->par[i]), &(ma->tair[i]), \
                          &(ma->tsoil[i]), &(ma->vpd[i]), &(ma->co2[i]), \
                          &(ma->ndep[i]), &(ma->wind[i]), \
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

void read_subdaily_met_data_binary(char **argv, control *c, met_arrays *ma)
{
    FILE  *fp;
    float *data = NULL;
    char   line[STRING_LENGTH];
    int    i = 0, j = 0;
    int    nvars = 12;
    int    skipped_lines = 0;
    double current_yr, temp_HOD;
    long   file_len, cnt;
    double hod_dummy;

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

    /* work out how big the file is */
    file_len = c->nrows*c->ncols;

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

    if ((ma->wind = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind array\n");
		exit(EXIT_FAILURE);
    }

    if ((ma->press = (double *)calloc(file_len, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for press array\n");
		exit(EXIT_FAILURE);
    }


    /* Divide up met chunk */
    cnt = 0;
    c->num_years = 0;
    current_yr = data[0];
    for (i = 0; i < c->nrows * c->ncols; i += c->ncols) {

            ma->year[cnt] = data[i];
            ma->doy[cnt] = data[i+1];
            temp_HOD = data[i+2];
            ma->rain[cnt] = data[i+3];
            ma->par[cnt] = data[i+4];
            ma->tair[cnt] = data[i+5];
            ma->tsoil[cnt] = data[i+6];
            ma->vpd[cnt] = data[i+7];
            ma->co2[cnt] = data[i+8];
            ma->ndep[cnt] = data[i+9];
            ma->wind[cnt] = data[i+10];
            ma->press[cnt] = data[i+11];

            /* Build an array of the unique years as we loop over the input file */
            if (current_yr != ma->year[cnt]) {
                c->num_years++;
                current_yr = ma->year[cnt];
            }

            cnt++;
    }

    fclose(fp);
    return;
}
