#include "read_met_file.h"

void read_met_data(char **argv, control *c, met *m)
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
    while (fgets(line, STRING_LENGTH, fp) != NULL) {
        /* ignore comment line */
        if (*line == '#')
            continue;
        file_len++;
    }
    rewind(fp);
    c->num_days = file_len;

    /* allocate memory for meteorological arrays */
    if ((m->year = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for year array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->prjday = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for prjday array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->sw_rad = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for sw_rad array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->tair = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tair array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->rain = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for rain array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->tsoil = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tsoil array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->tam = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tam array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->tpm = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for tpm array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->vpd_am = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->vpd_pm = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd_pm array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->vpd_avg = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for vpd_avg array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->co2 = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for co2 array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->ndep = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for ndep array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->wind = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->press = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for press array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->par = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for par array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->wind_am = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->wind_pm = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for wind_pm array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->sw_rad_am = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for sw_rad_am array\n");
		exit(EXIT_FAILURE);
    }

    if ((m->sw_rad_pm = (double *)calloc(c->num_days, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for sw_rad_pm array\n");
		exit(EXIT_FAILURE);
    }

    current_yr = m->year[0];

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
                          %lf,%lf,%lf,\
                          %lf,%lf,%lf,\
                          %lf,%lf", \
                          &(m->year[i]), &(m->prjday[i]), &(m->sw_rad[i]), \
                          &(m->tair[i]), &(m->rain[i]), &(m->tsoil[i]), \
                          &(m->tam[i]), &(m->tpm[i]), &(m->vpd_am[i]), \
                          &(m->vpd_pm[i]), &(m->vpd_avg[i]), &(m->co2[i]), \
                          &(m->ndep[i]), &(m->wind[i]), &(m->press[i]), \
                          &(m->par[i]), &(m->wind_am[i]), &(m->wind_pm[i]), \
                          &(m->sw_rad_am[i]), &(m->sw_rad_pm[i])) != nvars) {
            fprintf(stderr, "%s: badly formatted input in met file on line %d %d\n", \
                    *argv, (int)i+1+skipped_lines, nvars);
            exit(EXIT_FAILURE);
        }

        /* Build an array of the unique years as we loop over the input file */
        if (current_yr != m->year[i]) {
            c->num_years++;
            current_yr = m->year[i];
        }
        i++;
    }
    fclose(fp);
    return;
}
