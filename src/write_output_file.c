/* ============================================================================
* Print output file (ascii/binary)
*
*
*
* NOTES:
*   Currently I have only implemented the ascii version
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   25.02.2015
*
* =========================================================================== */
#include "write_output_file.h"


void open_output_file(control *c, char *fname, FILE **fp) {
    *fp = fopen(fname, "w");
    if (*fp == NULL)
        prog_error("Error opening output file for write on line", __LINE__);
}


void write_output_header(control *c, FILE **fp) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    int ncols = c->ovars; /* this is hardwired, but obv should match below! */
    int nrows = c->total_num_days;

    /* Git version */
    fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "YEAR,DOY,");

    /* plant */
    fprintf(*fp, "CF,LAI,CB,CW,CR,");

    /* water*/
    fprintf(*fp, "BETA,SWC,TRANS,SOIL_EVAP,CAN_EVAP,RUNOFF,");

    /* C fluxes */
    fprintf(*fp, "NPP\n");

    if (c->output_ascii == FALSE) {
        fprintf(*fp, "nrows=%d\n", nrows);
        fprintf(*fp, "ncols=%d\n", ncols);
    }

    return;
}

void write_daily_outputs_ascii(control *c, fluxes *f, state *s, int year,
                               int doy) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    float tonnes_per_ha_to_g_m2 = 100.0;
    int half_yr, offset;

    /* time stuff */
    fprintf(c->ofp, "%.10f,%.10f,", (double)year, (double)doy);


    /* plant */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f",
                    s->shoot, s->lai, s->stem,s->branch, s->root);

    /*
    ** FLUXES
    */

    /* water */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
            s->wtfac_root,s->pawater_root,f->transpiration,f->soil_evap,
            f->canopy_evap,f->runoff);


    /* C fluxes */
    fprintf(c->ofp, "%.10f\n", f->npp);

    return;
}


void save_daily_outputs_binary(control *c, fluxes *f, state *s, int year,
                                int doy, double *odata, long ocnt) {

    /* save everything and do a single big dump at the end */
    odata[ocnt] = (double)year;
    odata[ocnt+1] = (double)doy;
    odata[ocnt+2] = s->shoot;
    odata[ocnt+3] = s->lai;
    odata[ocnt+4] = s->branch;
    odata[ocnt+5] = s->stem;
    odata[ocnt+6] = s->root;
    odata[ocnt+7] = s->wtfac_root;
    odata[ocnt+8] = s->pawater_root;
    odata[ocnt+9] = f->transpiration;
    odata[ocnt+10] = f->soil_evap;
    odata[ocnt+11] = f->canopy_evap;
    odata[ocnt+12] = f->runoff;
    odata[ocnt+13] = f->npp;

    return;
}


int write_final_state(control *c, params *p, state *s)
{
    /*
    Write the final state to the input param file so we can easily restart
    the model. This function copies the input param file with the exception
    of anything in the git hash and the state which it replaces with the updated
    stuff.

    */

    char line[STRING_LENGTH];
    char saved_line[STRING_LENGTH];
    char section[STRING_LENGTH] = "";
    char prev_name[STRING_LENGTH] = "";
    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;
    int line_number = 0;
    int match = FALSE;

    while (fgets(line, sizeof(line), c->ifp) != NULL) {
        strcpy(saved_line, line);
        line_number++;
        start = lskip(rstrip(line));
        if (*start == ';' || *start == '#') {
            /* Per Python ConfigParser, allow '#' comments at start of line */
        }
        else if (*start == '[') {
            /* A "[section]" line */
            end = find_char_or_comment(start + 1, ']');
            if (*end == ']') {
                *end = '\0';
                strncpy0(section, start + 1, sizeof(section));
                *prev_name = '\0';

            }
            else if (!error) {
                /* No ']' found on section line */
                error = line_number;

            }
        }
        else if (*start && *start != ';') {
            /* Not a comment, must be a name[=:]value pair */
            end = find_char_or_comment(start, '=');
            if (*end != '=') {
                end = find_char_or_comment(start, ':');
            }
            if (*end == '=' || *end == ':') {
                *end = '\0';
                name = rstrip(start);
                value = lskip(end + 1);
                end = find_char_or_comment(value, '\0');
                if (*end == ';')
                    *end = '\0';
                rstrip(value);

                /* Valid name[=:]value pair found, call handler */
                strncpy0(prev_name, name, sizeof(prev_name));

                if (!ohandler(section, name, value, c, p, s, &match) && !error)
                    error = line_number;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = line_number;
                break;
            }
        }
        if (match == FALSE)
            fprintf(c->ofp, "%s", saved_line);
        else
            match = FALSE; /* reset match flag */
    }
    return error;

}


int ohandler(char *section, char *name, char *value, control *c, params *p,
             state *s, int *match)
{
    /*
    Search for matches of the git and state values and where found write the
    current state values to the output parameter file.

    - also added previous ncd as this potential can be changed internally
    */

    #define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    /*
    ** GIT
    */
    if (MATCH("git", "git_hash")) {
        fprintf(c->ofp, "git_hash = %s\n", c->git_code_ver);
        *match = TRUE;
    }

    /*
    ** PARAMS
    */
    if (MATCH("params", "previous_ncd")) {
        fprintf(c->ofp, "previous_ncd = %.10f\n", p->previous_ncd);
        *match = TRUE;
    }

    /*
    ** STATE
    */

    if (MATCH("state", "activesoil")) {
        fprintf(c->ofp, "activesoil = %.10f\n", s->activesoil);
        *match = TRUE;
    } else if (MATCH("state", "activesoiln")) {
        fprintf(c->ofp, "activesoiln = %.10f\n", s->activesoiln);
        *match = TRUE;
    } else if (MATCH("state", "age")) {
        fprintf(c->ofp, "age = %.10f\n", s->age);
        *match = TRUE;
    } else if (MATCH("state", "avg_albranch")) {
        fprintf(c->ofp, "avg_albranch = %.10f\n", s->avg_albranch);
        *match = TRUE;
    } else if (MATCH("state", "avg_alcroot")) {
        fprintf(c->ofp, "avg_alcroot = %.10f\n", s->avg_alcroot);
        *match = TRUE;
    } else if (MATCH("state", "avg_alleaf")) {
        fprintf(c->ofp, "avg_alleaf = %.10f\n", s->avg_alleaf);
        *match = TRUE;
    } else if (MATCH("state", "avg_alroot")) {
        fprintf(c->ofp, "avg_alroot = %.10f\n", s->avg_alroot);
        *match = TRUE;
    } else if (MATCH("state", "avg_alstem")) {
        fprintf(c->ofp, "avg_alstem = %.10f\n", s->avg_alstem);
        *match = TRUE;
    } else if (MATCH("state", "branch")) {
        fprintf(c->ofp, "branch = %.10f\n", s->branch);
        *match = TRUE;
    } else if (MATCH("state", "branchn")) {
        fprintf(c->ofp, "branchn = %.10f\n", s->branchn);
        *match = TRUE;
    } else if (MATCH("state", "canht")) {
        fprintf(c->ofp, "canht = %.10f\n", s->canht);
        *match = TRUE;
    } else if (MATCH("state", "croot")) {
        fprintf(c->ofp, "croot = %.10f\n", s->croot);
        *match = TRUE;
    } else if (MATCH("state", "crootn")) {
        fprintf(c->ofp, "crootn = %.10f\n", s->crootn);
        *match = TRUE;
    } else if (MATCH("state", "cstore")) {
        fprintf(c->ofp, "cstore = %.10f\n", s->cstore);
        *match = TRUE;
    } else if (MATCH("state", "inorgn")) {
        fprintf(c->ofp, "inorgn = %.10f\n", s->inorgn);
        *match = TRUE;
    } else if (MATCH("state", "lai")) {
        fprintf(c->ofp, "lai = %.10f\n", s->lai);
        *match = TRUE;
    } else if (MATCH("state", "metabsoil")) {
        fprintf(c->ofp, "metabsoil = %.10f\n", s->metabsoil);
        *match = TRUE;
    } else if (MATCH("state", "metabsoiln")) {
        fprintf(c->ofp, "metabsoiln = %.10f\n", s->metabsoiln);
        *match = TRUE;
    } else if (MATCH("state", "metabsurf")) {
        fprintf(c->ofp, "metabsurf = %.10f\n", s->metabsurf);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfn")) {
        fprintf(c->ofp, "metabsurfn = %.10f\n", s->metabsurfn);
        *match = TRUE;
    } else if (MATCH("state", "nstore")) {
        fprintf(c->ofp, "nstore = %.10f\n", s->nstore);
        *match = TRUE;
    } else if (MATCH("state", "passivesoil")) {
        fprintf(c->ofp, "passivesoil = %.10f\n", s->passivesoil);
        *match = TRUE;
    } else if (MATCH("state", "passivesoiln")) {
        fprintf(c->ofp, "passivesoiln = %.10f\n", s->passivesoiln);
        *match = TRUE;
    } else if (MATCH("state", "pawater_root")) {
        fprintf(c->ofp, "pawater_root = %.10f\n", s->pawater_root);
        *match = TRUE;
    } else if (MATCH("state", "pawater_topsoil")) {
        fprintf(c->ofp, "pawater_topsoil = %.10f\n", s->pawater_topsoil);
        *match = TRUE;
    } else if (MATCH("state", "prev_sma")) {
        fprintf(c->ofp, "prev_sma = %.10f\n", s->prev_sma);
        *match = TRUE;
    } else if (MATCH("state", "root")) {
        fprintf(c->ofp, "root = %.10f\n", s->root);
        *match = TRUE;
    } else if (MATCH("state", "root_depth")) {
        fprintf(c->ofp, "root_depth = %.10f\n", s->root_depth);
        *match = TRUE;
    } else if (MATCH("state", "rootn")) {
        fprintf(c->ofp, "rootn = %.10f\n", s->rootn);
        *match = TRUE;
    } else if (MATCH("state", "sapwood")) {
        fprintf(c->ofp, "sapwood = %.10f\n", s->sapwood);
        *match = TRUE;
    } else if (MATCH("state", "shoot")) {
        fprintf(c->ofp, "shoot = %.10f\n", s->shoot);
        *match = TRUE;
    } else if (MATCH("state", "shootn")) {
        fprintf(c->ofp, "shootn = %.10f\n", s->shootn);
        *match = TRUE;
    } else if (MATCH("state", "sla")) {
        fprintf(c->ofp, "sla = %.10f\n", s->sla);
        *match = TRUE;
    } else if (MATCH("state", "slowsoil")) {
        fprintf(c->ofp, "slowsoil = %.10f\n", s->slowsoil);
        *match = TRUE;
    } else if (MATCH("state", "slowsoiln")) {
        fprintf(c->ofp, "slowsoiln = %.10f\n", s->slowsoiln);
        *match = TRUE;
    } else if (MATCH("state", "stem")) {
        fprintf(c->ofp, "stem = %.10f\n", s->stem);
        *match = TRUE;
    } else if (MATCH("state", "stemn")) {
        fprintf(c->ofp, "stemn = %.10f\n", s->stemn);
        *match = TRUE;
    } else if (MATCH("state", "stemnimm")) {
        fprintf(c->ofp, "stemnimm = %.10f\n", s->stemnimm);
        *match = TRUE;
    } else if (MATCH("state", "stemnmob")) {
        fprintf(c->ofp, "stemnmob = %.10f\n", s->stemnmob);
        *match = TRUE;
    } else if (MATCH("state", "structsoil")) {
        fprintf(c->ofp, "structsoil = %.10f\n", s->structsoil);
        *match = TRUE;
    } else if (MATCH("state", "structsoiln")) {
        fprintf(c->ofp, "structsoiln = %.10f\n", s->structsoiln);
        *match = TRUE;
    } else if (MATCH("state", "structsurf")) {
        fprintf(c->ofp, "structsurf = %.10f\n", s->structsurf);
        *match = TRUE;
    } else if (MATCH("state", "structsurfn")) {
        fprintf(c->ofp, "structsurfn = %.10f\n", s->structsurfn);
        *match = TRUE;
    }

    return (1);
}
