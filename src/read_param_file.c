/* Read the .ini file into the various structures. */

#include "read_param_file.h"

int parse_ini_file(control *c, params *p, state *s)
{
    /*

    Loop through the file which is passed on the standard in, and break
    the file up into the relevant sections...

    */

    char line[STRING_LENGTH];
    char section[STRING_LENGTH] = "";
    char prev_name[STRING_LENGTH] = "";
    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;
    int line_number = 0;

    if ((c->ifp = fopen(c->cfg_fname, "r")) == NULL){
        prog_error("Error opening output file for write on line", __LINE__);
    }

    while (fgets(line, sizeof(line), c->ifp) != NULL) {
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

                if (!handler(section, name, value, c, p, s) && !error)
                    error = line_number;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = line_number;
                break;
            }
        }
    }

    if (c->print_options == END) {
        /* we need to re-read this file to dump the final state */
        rewind(c->ifp);
    }

    return error;


}



int handler(char *section, char *name, char *value, control *c,
            params *p, state *s)
{
    /*

    Assigns the values from the .INI file straight into the various
    structures

    */
    char *temp = value;

    #define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    /*
    ** GIT
    */
    if (MATCH("git", "git_hash")) {
        strcpy(c->git_hash, temp);
    }

    /*
    ** FILES
    */
    if (MATCH("files", "cfg_fname")) {
        /* remove quotation marks around the string
	    temp++;  removes first quote
	    temp[strlen(temp)-1] = 0;  removes last quote */
        strcpy(c->cfg_fname, temp);
    } else if (MATCH("files", "met_fname")) {
        strcpy(c->met_fname, temp);
    } else if (MATCH("files", "out_fname")) {
        strcpy(c->out_fname, temp);
    } else if (MATCH("files", "out_fname_hdr")) {
        strcpy(c->out_fname_hdr, temp);
    } else if (MATCH("files", "out_param_fname")) {
        strcpy(c->out_param_fname, temp);
    }

    /*
    ** CONTROL
    */
    if (MATCH("control", "alloc_model")) {
        if (strcmp(temp, "FIXED") == 0||
            strcmp(temp, "fixed") == 0)
            c->alloc_model = FIXED;
        else if (strcmp(temp, "GRASSES") == 0||
                 strcmp(temp, "grasses") == 0)
            c->alloc_model = GRASSES;
        else if (strcmp(temp, "ALLOMETRIC") == 0 ||
                 strcmp(temp, "allometric") == 0)
            c->alloc_model = ALLOMETRIC;
        else {
            fprintf(stderr, "Unknown alloc model: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "assim_model")) {
        if (strcmp(temp, "BEWDY") == 0||
            strcmp(temp, "bewdy") == 0)
            c->assim_model = BEWDY;
        else if (strcmp(temp, "MATE") == 0||
                 strcmp(temp, "mate") == 0)
            c->assim_model = MATE;
        else {
            fprintf(stderr, "Unknown photosynthesis model: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "calc_sw_params")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->calc_sw_params = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->calc_sw_params = TRUE;
        else {
            fprintf(stderr, "Unknown SW param option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "deciduous_model")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->deciduous_model = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->deciduous_model = TRUE;
        else {
            fprintf(stderr, "Unknown deciduous option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "fixed_stem_nc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixed_stem_nc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->fixed_stem_nc = TRUE;
        else {
            fprintf(stderr, "Unknown fixed_stem_nc option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "fixleafnc")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->fixleafnc = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->fixleafnc = TRUE;
        else {
            fprintf(stderr, "Unknown fixleafnc option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "grazing")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->grazing = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->grazing = TRUE;
        else {
            fprintf(stderr, "Unknown grazing option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "gs_model")) {
        if (strcmp(temp, "MEDLYN") == 0||
            strcmp(temp, "medlyn") == 0)
            c->gs_model = MEDLYN;
        else {
            fprintf(stderr, "Unknown gs model: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "model_optroot")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->model_optroot = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->model_optroot = TRUE;
        else {
            fprintf(stderr, "Unknown model_optroot option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "modeljm")) {
        c->modeljm = atoi(value);
    } else if (MATCH("control", "ncycle")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->ncycle = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->ncycle = TRUE;
        else {
            fprintf(stderr, "Unknown ncycle option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "nuptake_model")) {
        c->nuptake_model = atoi(value);
    } else if (MATCH("control", "output_ascii")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->output_ascii = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->output_ascii = TRUE;
        else {
            fprintf(stderr, "Unknown output_ascii option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "passiveconst")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->passiveconst = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->passiveconst = TRUE;
        else {
            fprintf(stderr, "Unknown passiveconst option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "print_options")) {
        if (strcmp(temp, "Daily") == 0 ||
            strcmp(temp, "DAILY") == 0 ||
            strcmp(temp, "daily") == 0)
            c->print_options = DAILY;
        else if (strcmp(temp, "End") == 0 ||
            strcmp(temp, "END") == 0 ||
            strcmp(temp, "end") == 0)
            c->print_options = END;
        else {
            fprintf(stderr, "Unknown print option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    } else if (MATCH("control", "ps_pathway")) {
        if (strcmp(temp, "C3") == 0 ||
            strcmp(temp, "c3") == 0)
            c->ps_pathway = C3;
        else if (strcmp(temp, "End") == 0 ||
            strcmp(temp, "C4") == 0 ||
            strcmp(temp, "c4") == 0)
            c->ps_pathway = C4;
        else {
            fprintf(stderr, "Unknown ps pathway : %s\n", temp);
            exit(EXIT_FAILURE);
        }
     } else if (MATCH("control", "respiration_model")) {
         if (strcmp(temp, "FIXED") == 0||
             strcmp(temp, "fixed") == 0)
             c->respiration_model = FIXED;
         else if (strcmp(temp, "TEMPERATURE") == 0||
             strcmp(temp, "temperature") == 0)
             c->respiration_model = TEMPERATURE;
         else if (strcmp(temp, "BIOMASS") == 0||
             strcmp(temp, "biomass") == 0)
                 c->respiration_model = BIOMASS;
         else {
             fprintf(stderr, "Unknown respiration model: %s\n", temp);
             exit(EXIT_FAILURE);
         }
    } else if (MATCH("control", "strfloat")) {
        c->strfloat = atoi(value);
        /*if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0)
            c->strfloat = FALSE;
        else if (strcmp(temp, "True") == 0 ||
            strcmp(temp, "TRUE") == 0 ||
            strcmp(temp, "true") == 0)
            c->strfloat = TRUE;
        else {
            fprintf(stderr, "Unknown strfloat option: %s\n", temp);
            exit(EXIT_FAILURE);
        }*/
    } else if (MATCH("control", "sw_stress_model")) {
        c->sw_stress_model = atoi(value);
    } else if (MATCH("control", "trans_model")) {
        c->trans_model = atoi(value);
    } else if (MATCH("control", "use_eff_nc")) {
        c->use_eff_nc = atoi(value);
    } else if (MATCH("control", "water_stress")) {
        if (strcmp(temp, "False") == 0 ||
            strcmp(temp, "FALSE") == 0 ||
            strcmp(temp, "false") == 0) {
            c->water_stress = FALSE;
            fprintf(stderr, "\nYou have turned off the drought stress??\n");
        } else if (strcmp(temp, "True") == 0 ||
                   strcmp(temp, "TRUE") == 0 ||
                   strcmp(temp, "true") == 0) {
            c->water_stress = TRUE;
        } else {
            fprintf(stderr, "Unknown water stress option: %s\n", temp);
            exit(EXIT_FAILURE);
        }
    }



    /*
    ** State
    */
    if (MATCH("state", "activesoil")) {
        s->activesoil = atof(value);
    } else if (MATCH("state", "activesoiln")) {
        s->activesoiln = atof(value);
    } else if (MATCH("state", "age")) {
        s->age = atof(value);
    } else if (MATCH("state", "avg_albranch")) {
        s->avg_albranch = atof(value);
    } else if (MATCH("state", "avg_alcroot")) {
        s->avg_alcroot = atof(value);
    } else if (MATCH("state", "avg_alleaf")) {
        s->avg_alleaf = atof(value);
    } else if (MATCH("state", "avg_alroot")) {
        s->avg_alroot = atof(value);
    } else if (MATCH("state", "avg_alstem")) {
        s->avg_alstem = atof(value);
    } else if (MATCH("state", "branch")) {
        s->branch = atof(value);
    } else if (MATCH("state", "branchn")) {
        s->branchn = atof(value);
    } else if (MATCH("state", "canht")) {
        s->canht = atof(value);
    } else if (MATCH("state", "croot")) {
        s->croot = atof(value);
    } else if (MATCH("state", "crootn")) {
        s->crootn = atof(value);
    } else if (MATCH("state", "cstore")) {
        s->cstore = atof(value);
    } else if (MATCH("state", "inorgn")) {
        s->inorgn = atof(value);
    } else if (MATCH("state", "lai")) {
        s->lai = atof(value);
    } else if (MATCH("state", "max_lai")) {
        s->max_lai = atof(value);
    } else if (MATCH("state", "max_shoot")) {
        s->max_shoot = atof(value);
    } else if (MATCH("state", "metabsoil")) {
        s->metabsoil = atof(value);
    } else if (MATCH("state", "metabsoiln")) {
        s->metabsoiln = atof(value);
    } else if (MATCH("state", "metabsurf")) {
        s->metabsurf = atof(value);
    } else if (MATCH("state", "metabsurfn")) {
        s->metabsurfn = atof(value);
    } else if (MATCH("state", "nstore")) {
        s->nstore = atof(value);
    } else if (MATCH("state", "passivesoil")) {
        s->passivesoil = atof(value);
    } else if (MATCH("state", "passivesoiln")) {
        s->passivesoiln = atof(value);
    } else if (MATCH("state", "pawater_root")) {
        s->pawater_root = atof(value);
    } else if (MATCH("state", "pawater_topsoil")) {
        s->pawater_topsoil = atof(value);
    } else if (MATCH("state", "prev_sma")) {
        s->prev_sma = atof(value);
    } else if (MATCH("state", "root")) {
        s->root = atof(value);
    } else if (MATCH("state", "root_depth")) {
        s->root_depth = atof(value);
    } else if (MATCH("state", "rootn")) {
        s->rootn = atof(value);
    } else if (MATCH("state", "sapwood")) {
        s->sapwood = atof(value);
    } else if (MATCH("state", "shoot")) {
        s->shoot = atof(value);
    } else if (MATCH("state", "shootn")) {
        s->shootn = atof(value);
    } else if (MATCH("state", "sla")) {
        s->sla = atof(value);
    } else if (MATCH("state", "slowsoil")) {
        s->slowsoil = atof(value);
    } else if (MATCH("state", "slowsoiln")) {
        s->slowsoiln = atof(value);
    } else if (MATCH("state", "stem")) {
        s->stem = atof(value);
    } else if (MATCH("state", "stemn")) {
        s->stemn = atof(value);
    } else if (MATCH("state", "stemnimm")) {
        s->stemnimm = atof(value);
    } else if (MATCH("state", "stemnmob")) {
        s->stemnmob = atof(value);
    } else if (MATCH("state", "structsoil")) {
        s->structsoil = atof(value);
    } else if (MATCH("state", "structsoiln")) {
        s->structsoiln = atof(value);
    } else if (MATCH("state", "structsurf")) {
        s->structsurf = atof(value);
    } else if (MATCH("state", "structsurfn")) {
        s->structsurfn = atof(value);
    }

    /* Params */
    if (MATCH("params", "actncmax")) {
        p->actncmax = atof(value);
    } else if (MATCH("params", "actncmin")) {
        p->actncmin = atof(value);
    } else if (MATCH("params", "adapt")) {
        p->adapt = atof(value);
    } else if (MATCH("params", "ageold")) {
        p->ageold = atof(value);
    } else if (MATCH("params", "ageyoung")) {
        p->ageyoung = atof(value);
    } else if (MATCH("params", "albedo")) {
        p->albedo = atof(value);
    } else if (MATCH("params", "alpha_c4")) {
        p->alpha_c4 = atof(value);
    } else if (MATCH("params", "alpha_j")) {
        p->alpha_j = atof(value);
    } else if (MATCH("params", "b_root")) {
        p->b_root = atof(value);
    } else if (MATCH("params", "b_topsoil")) {
        p->b_topsoil = atof(value);
    } else if (MATCH("params", "bdecay")) {
        p->bdecay = atof(value);
    } else if (MATCH("params", "branch0")) {
        p->branch0 = atof(value);
    } else if (MATCH("params", "branch1")) {
        p->branch1 = atof(value);
    } else if (MATCH("params", "c_alloc_bmax")) {
        p->c_alloc_bmax = atof(value);
    } else if (MATCH("params", "c_alloc_bmin")) {
        p->c_alloc_bmin = atof(value);
    } else if (MATCH("params", "c_alloc_cmax")) {
        p->c_alloc_cmax = atof(value);
    } else if (MATCH("params", "c_alloc_fmax")) {
        p->c_alloc_fmax = atof(value);
    } else if (MATCH("params", "c_alloc_fmin")) {
        p->c_alloc_fmin = atof(value);
    } else if (MATCH("params", "c_alloc_rmax")) {
        p->c_alloc_rmax = atof(value);
    } else if (MATCH("params", "c_alloc_rmin")) {
        p->c_alloc_rmin = atof(value);
    } else if (MATCH("params", "cfracts")) {
        p->cfracts = atof(value);
    } else if (MATCH("params", "crdecay")) {
        p->crdecay = atof(value);
    } else if (MATCH("params", "cretrans")) {
        p->cretrans = atof(value);
    } else if (MATCH("params", "croot0")) {
        p->croot0 = atof(value);
    } else if (MATCH("params", "croot1")) {
        p->croot1 = atof(value);
    } else if (MATCH("params", "ctheta_root")) {
        p->ctheta_root = atof(value);
    } else if (MATCH("params", "ctheta_topsoil")) {
        p->ctheta_topsoil = atof(value);
    } else if (MATCH("params", "cue")) {
        p->cue = atof(value);
    } else if (MATCH("params", "d0")) {
        p->d0 = atof(value);
    } else if (MATCH("params", "d0x")) {
        p->d0x = atof(value);
    } else if (MATCH("params", "d1")) {
        p->d1 = atof(value);
    } else if (MATCH("params", "delsj")) {
        p->delsj = atof(value);
    } else if (MATCH("params", "density")) {
        p->density = atof(value);
    } else if (MATCH("params", "direct_frac")) {
        p->direct_frac = atof(value);
    } else if (MATCH("params", "displace_ratio")) {
        p->displace_ratio = atof(value);
    } else if (MATCH("params", "disturbance_doy")) {
        p->disturbance_doy = atof(value);
    } else if (MATCH("params", "dz0v_dh")) {
        p->dz0v_dh = atof(value);
    } else if (MATCH("params", "eac")) {
        p->eac = atof(value);
    } else if (MATCH("params", "eag")) {
        p->eag = atof(value);
    } else if (MATCH("params", "eaj")) {
        p->eaj = atof(value);
    } else if (MATCH("params", "eao")) {
        p->eao = atof(value);
    } else if (MATCH("params", "eav")) {
        p->eav = atof(value);
    } else if (MATCH("params", "edj")) {
        p->edj = atof(value);
    } else if (MATCH("params", "faecescn")) {
        p->faecescn = atof(value);
    } else if (MATCH("params", "faecesn")) {
        p->faecesn = atof(value);
    } else if (MATCH("params", "fdecay")) {
        p->fdecay = atof(value);
    } else if (MATCH("params", "fdecaydry")) {
        p->fdecaydry = atof(value);
    } else if (MATCH("params", "fhw")) {
        p->fhw = atof(value);
    } else if (MATCH("params", "finesoil")) {
        p->finesoil = atof(value);
    } else if (MATCH("params", "fracfaeces")) {
        p->fracfaeces = atof(value);
    } else if (MATCH("params", "fracteaten")) {
        p->fracteaten = atof(value);
    } else if (MATCH("params", "fractosoil")) {
        p->fractosoil = atof(value);
    } else if (MATCH("params", "fractup_soil")) {
        p->fractup_soil = atof(value);
    } else if (MATCH("params", "fretrans")) {
        p->fretrans = atof(value);
    } else if (MATCH("params", "g1")) {
        p->g1 = atof(value);
    } else if (MATCH("params", "gamstar25")) {
        p->gamstar25 = atof(value);
    } else if (MATCH("params", "growth_efficiency")) {
        p->growth_efficiency = atof(value);
    } else if (MATCH("params", "gamstar25")) {
        p->gamstar25 = atof(value);
    } else if (MATCH("params", "height0")) {
        p->height0 = atof(value);
    } else if (MATCH("params", "height1")) {
        p->height1 = atof(value);
    } else if (MATCH("params", "heighto")) {
        p->heighto = atof(value);
    } else if (MATCH("params", "htpower")) {
        p->htpower = atof(value);
    } else if (MATCH("params", "intercep_frac")) {
        p->intercep_frac = atof(value);
    } else if (MATCH("params", "jmax")) {
        p->jmax = atof(value);
    } else if (MATCH("params", "jmaxna")) {
        p->jmaxna = atof(value);
    } else if (MATCH("params", "jmaxnb")) {
        p->jmaxnb = atof(value);
    } else if (MATCH("params", "jv_intercept")) {
        p->jv_intercept = atof(value);
    } else if (MATCH("params", "jv_slope")) {
        p->jv_slope = atof(value);
    } else if (MATCH("params", "kc25")) {
        p->kc25 = atof(value);
    } else if (MATCH("params", "kdec1")) {
        p->kdec1 = atof(value);
    } else if (MATCH("params", "kdec2")) {
        p->kdec2 = atof(value);
    } else if (MATCH("params", "kdec3")) {
        p->kdec3 = atof(value);
    } else if (MATCH("params", "kdec4")) {
        p->kdec4 = atof(value);
    } else if (MATCH("params", "kdec5")) {
        p->kdec5 = atof(value);
    } else if (MATCH("params", "kdec6")) {
        p->kdec6 = atof(value);
    } else if (MATCH("params", "kdec7")) {
        p->kdec7 = atof(value);
    } else if (MATCH("params", "knl")) {
        p->knl = atof(value);
    } else if (MATCH("params", "ko25")) {
        p->ko25 = atof(value);
    } else if (MATCH("params", "kq10")) {
        p->kq10 = atof(value);
    } else if (MATCH("params", "kr")) {
        p->kr = atof(value);
    } else if (MATCH("params", "lai_closed")) {
        p->lai_closed = atof(value);
    } else if (MATCH("params", "latitude")) {
        p->latitude = atof(value);
    } else if (MATCH("params", "leafsap0")) {
        p->leafsap0 = atof(value);
    } else if (MATCH("params", "leafsap1")) {
        p->leafsap1 = atof(value);
    } else if (MATCH("params", "ligfaeces")) {
        p->ligfaeces = atof(value);
    } else if (MATCH("params", "ligroot")) {
        p->ligroot = atof(value);
    } else if (MATCH("params", "ligshoot")) {
        p->ligshoot = atof(value);
    } else if (MATCH("params", "liteffnc")) {
        p->liteffnc = atof(value);
    } else if (MATCH("params", "max_intercep_lai")) {
        p->max_intercep_lai = atof(value);
    } else if (MATCH("params", "measurement_temp")) {
        p->measurement_temp = atof(value);
    } else if (MATCH("params", "ncbnew")) {
        p->ncbnew = atof(value);
    } else if (MATCH("params", "ncbnewz")) {
        p->ncbnewz = atof(value);
    } else if (MATCH("params", "nccnew")) {
        p->nccnew = atof(value);
    } else if (MATCH("params", "nccnewz")) {
        p->nccnewz = atof(value);
    } else if (MATCH("params", "ncmaxfold")) {
        p->ncmaxfold = atof(value);
    } else if (MATCH("params", "ncmaxfyoung")) {
        p->ncmaxfyoung = atof(value);
    } else if (MATCH("params", "ncmaxr")) {
        p->ncmaxr = atof(value);
    } else if (MATCH("params", "ncrfac")) {
        p->ncrfac = atof(value);
    } else if (MATCH("params", "ncwimm")) {
        p->ncwimm = atof(value);
    } else if (MATCH("params", "ncwimmz")) {
        p->ncwimmz = atof(value);
    } else if (MATCH("params", "ncwnew")) {
        p->ncwnew = atof(value);
    } else if (MATCH("params", "ncwnewz")) {
        p->ncwnewz = atof(value);
    } else if (MATCH("params", "nf_crit")) {
        p->nf_crit = atof(value);
    } else if (MATCH("params", "nf_min")) {
        p->nf_min = atof(value);
    } else if (MATCH("params", "nmax")) {
        p->nmax = atof(value);
    } else if (MATCH("params", "nmin")) {
        p->nmin = atof(value);
    } else if (MATCH("params", "nmin0")) {
        p->nmin0 = atof(value);
    } else if (MATCH("params", "nmincrit")) {
        p->nmincrit = atof(value);
    } else if (MATCH("params", "ntheta_root")) {
        p->ntheta_root = atof(value);
    } else if (MATCH("params", "ntheta_topsoil")) {
        p->ntheta_topsoil = atof(value);
    } else if (MATCH("params", "nuptakez")) {
        p->nuptakez = atof(value);
    } else if (MATCH("params", "oi")) {
        p->oi = atof(value);
    } else if (MATCH("params", "passivesoilnz")) {
        p->passivesoilnz = atof(value);
    } else if (MATCH("params", "passivesoilz")) {
        p->passivesoilz = atof(value);
    } else if (MATCH("params", "passncmax")) {
        p->passncmax = atof(value);
    } else if (MATCH("params", "passncmin")) {
        p->passncmin = atof(value);
    } else if (MATCH("params", "prescribed_leaf_NC")) {
        p->prescribed_leaf_NC = atof(value);
    } else if (MATCH("params", "previous_ncd")) {
        p->previous_ncd = atof(value);
    } else if (MATCH("params", "psi_sat_root")) {
        p->psi_sat_root = atof(value);
    } else if (MATCH("params", "psi_sat_topsoil")) {
        p->psi_sat_topsoil = atof(value);
    } else if (MATCH("params", "qs")) {
        p->qs = atof(value);
    } else if (MATCH("params", "r0")) {
        p->r0 = atof(value);
    } else if (MATCH("params", "rateloss")) {
        p->rateloss = atof(value);
    } else if (MATCH("params", "rateuptake")) {
        p->rateuptake = atof(value);
    } else if (MATCH("params", "rdecay")) {
        p->rdecay = atof(value);
    } else if (MATCH("params", "rdecaydry")) {
        p->rdecaydry = atof(value);
    } else if (MATCH("params", "retransmob")) {
        p->retransmob = atof(value);
    } else if (MATCH("params", "rfmult")) {
        p->rfmult = atof(value);
    } else if (MATCH("params", "rooting_depth")) {
        p->rooting_depth = atof(value);
    } else if (MATCH("params", "rootsoil_type")) {
        strcpy(p->rootsoil_type, value);
    } else if (MATCH("params", "rretrans")) {
        p->rretrans = atof(value);
    } else if (MATCH("params", "sapturnover")) {
        p->sapturnover = atof(value);
    } else if (MATCH("params", "sla")) {
        p->sla = atof(value);
    } else if (MATCH("params", "slamax")) {
        p->slamax = atof(value);
    } else if (MATCH("params", "slazero")) {
        p->slazero = atof(value);
    } else if (MATCH("params", "slowncmax")) {
        p->slowncmax = atof(value);
    } else if (MATCH("params", "slowncmin")) {
        p->slowncmin = atof(value);
    } else if (MATCH("params", "store_transfer_len")) {
        p->store_transfer_len = atof(value);
    } else if (MATCH("params", "structcn")) {
        p->structcn = atof(value);
    } else if (MATCH("params", "structrat")) {
        p->structrat = atof(value);
    } else if (MATCH("params", "targ_sens")) {
        p->targ_sens = atof(value);
    } else if (MATCH("params", "theta")) {
        p->theta = atof(value);
    } else if (MATCH("params", "theta_sat_root")) {
        p->theta_sat_root = atof(value);
    } else if (MATCH("params", "theta_sat_topsoil")) {
        p->theta_sat_topsoil = atof(value);
    } else if (MATCH("params", "topsoil_depth")) {
        p->topsoil_depth = atof(value);
    } else if (MATCH("params", "topsoil_type")) {
        strcpy(p->topsoil_type, value);
    } else if (MATCH("params", "vcmax")) {
        p->vcmax = atof(value);
    } else if (MATCH("params", "vcmaxna")) {
        p->vcmaxna = atof(value);
    } else if (MATCH("params", "vcmaxnb")) {
        p->vcmaxnb = atof(value);
    } else if (MATCH("params", "watdecaydry")) {
        p->watdecaydry = atof(value);
    } else if (MATCH("params", "watdecaywet")) {
        p->watdecaywet = atof(value);
    } else if (MATCH("params", "wcapac_root")) {
        p->wcapac_root = atof(value);
    } else if (MATCH("params", "wcapac_topsoil")) {
        p->wcapac_topsoil = atof(value);
    } else if (MATCH("params", "wdecay")) {
        p->wdecay = atof(value);
    } else if (MATCH("params", "wetloss")) {
        p->wetloss = atof(value);
    } else if (MATCH("params", "wretrans")) {
        p->wretrans = atof(value);
    } else if (MATCH("params", "z0h_z0m")) {
        p->z0h_z0m = atof(value);
    }

    return (1);
}
