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

void write_output_subdaily_header(control *c, FILE **fp) {
    /*
        Write 30 min fluxes headers to an output CSV file. This is very basic
        for now...
    */

    /* Git version */
    fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "year,doy,hod,");

    /*
    ** Canopy stuff...
    */
    fprintf(*fp, "an_canopy,rd_canopy,gsc_canopy,");
    fprintf(*fp, "apar_canopy,trans_canopy,tleaf\n");
    return;
}

void write_output_header(control *c, FILE **fp) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    int ncols = 86;  /* change with number of variables ? total count below is 93 */
    int nrows = c->num_days;

    /* Git version */
    fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "year,doy,");

    /*
    ** STATE
    */

    /* water*/
    fprintf(*fp, "wtfac_root,wtfac_topsoil,pawater_root,");

    /* plant */
    fprintf(*fp, "shoot,lai,branch,stem,root,croot,");
    fprintf(*fp, "shootn,branchn,stemn,rootn,crootn,");
    fprintf(*fp, "shootp,branchp,stemp,rootp,crootp,");
    fprintf(*fp, "cstore,nstore,pstore,");

    /* belowground */
    fprintf(*fp, "soilc,soiln,soilp,inorgn,");
    fprintf(*fp, "inorgp,inorgavlp,inorglabp,inorgsorbp,inorgssorbp,inorgoccp,inorgparp,");
    fprintf(*fp, "litterc,littercag,littercbg,litternag,litternbg,");
    fprintf(*fp, "litterpag,litterpbg,");
    fprintf(*fp, "activesoil,slowsoil,passivesoil,");
    fprintf(*fp, "activesoiln,slowsoiln,passivesoiln,activesoilp,slowsoilp,passivesoilp,");

    /*
    ** FLUXES
    */

    /* water */
    fprintf(*fp, "et,transpiration,soil_evap,canopy_evap,runoff,");
    fprintf(*fp, "gs_mol_m2_sec,ga_mol_m2_sec,");

    /* litter */
    fprintf(*fp, "deadleaves,deadbranch,deadstems,deadroots,deadcroots,");
    fprintf(*fp, "deadleafn,deadbranchn,deadstemn,deadrootn,deadcrootn,");
    fprintf(*fp, "deadleafp,deadbranchp,deadstemp,deadrootp,deadcrootp,");
    

    /* C fluxes */
    fprintf(*fp, "nep,gpp,npp,hetero_resp,auto_resp,apar,");

    /* C, N and P growth */
    fprintf(*fp, "cpleaf,cpbranch,cpstem,cproot,cpcroot,");
    fprintf(*fp, "npleaf,npbranch,npstemimm,npstemmob,nproot,npcroot,");
    fprintf(*fp, "ppleaf,ppbranch,ppstemimm,ppstemmob,pproot,ppcroot,");
    

    /* N stuff */
    fprintf(*fp, "nuptake,ngross,nmineralisation,nloss,");
    
    /* P stuff */
    fprintf(*fp, "puptake,pgross,pmineralisation,ploss,");

    /* traceability stuff */
    fprintf(*fp, "tfac_soil_decomp,c_into_active,c_into_slow,");
    fprintf(*fp, "c_into_passive,active_to_slow,active_to_passive,");
    fprintf(*fp, "slow_to_active,slow_to_passive,passive_to_active,");
    fprintf(*fp, "co2_rel_from_surf_struct_litter,");
    fprintf(*fp, "co2_rel_from_soil_struct_litter,");
    fprintf(*fp, "co2_rel_from_surf_metab_litter,");
    fprintf(*fp, "co2_rel_from_soil_metab_litter,");
    fprintf(*fp, "co2_rel_from_active_pool,");
    fprintf(*fp, "co2_rel_from_slow_pool,");
    fprintf(*fp, "co2_rel_from_passive_pool,");

    /* extra priming stuff */
    fprintf(*fp, "root_exc,");
    fprintf(*fp, "root_exn,");
    fprintf(*fp, "root_exp,");
    fprintf(*fp, "co2_released_exud,");
    fprintf(*fp, "factive,");
    fprintf(*fp, "rtslow,");
    fprintf(*fp, "rexc_cue,");

    /* Misc */
    fprintf(*fp, "leafretransn,"); 
    fprintf(*fp, "leafretransp\n");
    fprintf(*fp, "predawn_swp,");
    fprintf(*fp, "midday_lwp,");
    fprintf(*fp, "leafretransn\n");


    if (c->output_ascii == FALSE) {
        fprintf(*fp, "nrows=%d\n", nrows);
        fprintf(*fp, "ncols=%d\n", ncols);
    }
    return;
}

void write_subdaily_outputs_ascii(control *c, canopy_wk *cw, double year,
                                  double doy, int hod) {
    /*
        Write sub-daily canopy fluxes - very basic for now
    */

    /* time stuff */
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f,", year, doy, (double)hod);

    /* Canopy stuff */
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f,",
                       cw->an_canopy, cw->rd_canopy, cw->gsc_canopy);
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f\n",
                       cw->apar_canopy, cw->trans_canopy, cw->tleaf_new);

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


    /* time stuff */
    fprintf(c->ofp, "%.10f,%.10f,", (double)year, (double)doy);

    /*
    ** STATE

    */

    /* water*/
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
            s->wtfac_root,s->wtfac_topsoil,s->pawater_root);

    /* plant */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->shoot,s->lai,s->branch,s->stem,s->root,s->croot);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->shootn,s->branchn,s->stemn,s->rootn,s->crootn);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->shootp,s->branchp,s->stemp,s->rootp,s->crootp);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",   
                    s->cstore,s->nstore,s->pstore);

    /* belowground */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    s->soilc,s->soiln,s->soilp,s->inorgn);
                    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->inorgp,s->inorgavlp,s->inorglabp,s->inorgsorbp,s->inorgssorbp,
                    s->inorgoccp,s->inorgparp);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->litterc,s->littercag,s->littercbg,s->litternag,s->litternbg);
    
    fprintf(c->ofp, "%.10f,%.10f,",
                    s->litterpag,s->litterpbg);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    s->activesoil,s->slowsoil,s->passivesoil);
    
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->activesoiln,s->slowsoiln,s->passivesoiln,
                    s->activesoilp,s->slowsoilp,s->passivesoilp);
    /*
    ** FLUXES
    */

    /* water */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->et,f->transpiration,f->soil_evap,f->canopy_evap);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->runoff,f->gs_mol_m2_sec,f->ga_mol_m2_sec);

    /* litter */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->deadleaves,f->deadbranch,f->deadstems,f->deadroots,
                    f->deadcroots);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->deadleafn,f->deadbranchn,f->deadstemn,f->deadrootn,
                    f->deadcrootn);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->deadleafp,f->deadbranchp,f->deadstemp,f->deadrootp,
                    f->deadcrootp);

    /* C fluxes */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->nep,f->gpp,f->npp,f->hetero_resp,f->auto_resp,
                    f->apar);

    /* C N and P growth */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->cpleaf,f->cpbranch,f->cpstem,f->cproot,f->cpcroot);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->npleaf,f->npbranch,f->npstemimm,f->npstemmob,f->nproot,
                    f->npcroot);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->ppleaf,f->ppbranch,f->ppstemimm,f->ppstemmob,f->pproot,
                    f->ppcroot);

    /* N stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->nuptake,f->ngross,f->nmineralisation,f->nloss);
    
    /* P stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->puptake,f->pgross,f->pmineralisation,f->ploss);


    /* traceability stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->tfac_soil_decomp,f->c_into_active,f->c_into_slow);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->c_into_passive,f->active_to_slow,f->active_to_passive);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->slow_to_active,f->slow_to_passive,f->passive_to_active);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_active_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_slow_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_passive_pool);

    /* extra priming stuff */
    fprintf(c->ofp, "%.10f,", f->root_exc);
    fprintf(c->ofp, "%.10f,", f->root_exn);
    fprintf(c->ofp, "%.10f,", f->root_exp);
    fprintf(c->ofp, "%.10f,", f->co2_released_exud);
    fprintf(c->ofp, "%.10f,", f->factive);
    fprintf(c->ofp, "%.10f,", f->rtslow);
    fprintf(c->ofp, "%.10f,", f->rexc_cue);

    /* Misc */
    fprintf(c->ofp, "%.10f,", f->leafretransn);
    fprintf(c->ofp, "%.10f\n", f->leafretransp);
    fprintf(c->ofp, "%.10f,", s->predawn_swp);
    fprintf(c->ofp, "%.10f,", s->midday_lwp);
    fprintf(c->ofp, "%.10f\n", f->leafretransn);


    return;
}

void write_daily_outputs_binary(control *c, fluxes *f, state *s, int year,
                                int doy) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    double temp;

    /* time stuff */
    temp = (double)year;
    fwrite(&temp, sizeof(double), 1, c->ofp);
    temp = (double)doy;
    fwrite(&temp, sizeof(double), 1, c->ofp);


    /* plant */
    fwrite(&(s->shoot), sizeof(double), 1, c->ofp);
    fwrite(&(s->lai), sizeof(double), 1, c->ofp);
    fwrite(&(s->branch), sizeof(double), 1, c->ofp);
    fwrite(&(s->stem), sizeof(double), 1, c->ofp);
    fwrite(&(s->root), sizeof(double), 1, c->ofp);

    /* water */
    fwrite(&(s->wtfac_root), sizeof(double), 1, c->ofp);
    fwrite(&(s->pawater_root), sizeof(double), 1, c->ofp);
    fwrite(&(f->transpiration), sizeof(double), 1, c->ofp);
    fwrite(&(f->soil_evap), sizeof(double), 1, c->ofp);
    fwrite(&(f->canopy_evap), sizeof(double), 1, c->ofp);
    fwrite(&(f->runoff), sizeof(double), 1, c->ofp);

    /* C fluxes */
    fwrite(&(f->npp), sizeof(double), 1, c->ofp);

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
    } else if (MATCH("state", "activesoilp")) {
        fprintf(c->ofp, "activesoilp = %.10f\n", s->activesoilp);
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
    } else if (MATCH("state", "branchp")) {
        fprintf(c->ofp, "branchp = %.10f\n", s->branchp);
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
    } else if (MATCH("state", "crootp")) {
        fprintf(c->ofp, "crootp = %.10f\n", s->crootp);
        *match = TRUE;
    } else if (MATCH("state", "cstore")) {
        fprintf(c->ofp, "cstore = %.10f\n", s->cstore);
        *match = TRUE;
    } else if (MATCH("state", "inorgn")) {
        fprintf(c->ofp, "inorgn = %.10f\n", s->inorgn);
        *match = TRUE;
    } else if (MATCH("state", "inorgp")) {
        fprintf(c->ofp, "inorgp = %.10f\n", s->inorgp);
        *match = TRUE;
    } else if (MATCH("state", "inorgavlp")) {
        fprintf(c->ofp, "inorgavlp = %.10f\n", s->inorgavlp);
        *match = TRUE;
    } else if (MATCH("state", "inorglabp")) {
        fprintf(c->ofp, "inorglabp = %.10f\n", s->inorglabp);
        *match = TRUE;
    } else if (MATCH("state", "inorgsorbp")) {
        fprintf(c->ofp, "inorgsorbp = %.10f\n", s->inorgsorbp);
        *match = TRUE;
    } else if (MATCH("state", "inorgssorbp")) {
        fprintf(c->ofp, "inorgssorbp = %.10f\n", s->inorgssorbp);
        *match = TRUE;
    } else if (MATCH("state", "inorgoccp")) {
        fprintf(c->ofp, "inorgoccp = %.10f\n", s->inorgoccp);
        *match = TRUE;
    } else if (MATCH("state", "inorgparp")) {
        fprintf(c->ofp, "inorgparp = %.10f\n", s->inorgparp);
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
    } else if (MATCH("state", "metabsoilp")) {
        fprintf(c->ofp, "metabsoilp = %.10f\n", s->metabsoilp);
        *match = TRUE;
    } else if (MATCH("state", "metabsurf")) {
        fprintf(c->ofp, "metabsurf = %.10f\n", s->metabsurf);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfn")) {
        fprintf(c->ofp, "metabsurfn = %.10f\n", s->metabsurfn);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfp")) {
        fprintf(c->ofp, "metabsurfp = %.10f\n", s->metabsurfp);
        *match = TRUE;
    } else if (MATCH("state", "nstore")) {
        fprintf(c->ofp, "nstore = %.10f\n", s->nstore);
        *match = TRUE;
    } else if (MATCH("state", "pstore")) {
        fprintf(c->ofp, "pstore = %.10f\n", s->pstore);
        *match = TRUE;
    } else if (MATCH("state", "passivesoil")) {
        fprintf(c->ofp, "passivesoil = %.10f\n", s->passivesoil);
        *match = TRUE;
    } else if (MATCH("state", "passivesoiln")) {
        fprintf(c->ofp, "passivesoiln = %.10f\n", s->passivesoiln);
        *match = TRUE;
    } else if (MATCH("state", "passivesoilp")) {
        fprintf(c->ofp, "passivesoilp = %.10f\n", s->passivesoilp);
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
    } else if (MATCH("state", "rootp")) {
        fprintf(c->ofp, "rootp = %.10f\n", s->rootp);
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
    } else if (MATCH("state", "shootp")) {
        fprintf(c->ofp, "shootp = %.10f\n", s->shootp);
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
    } else if (MATCH("state", "slowsoilp")) {
        fprintf(c->ofp, "slowsoilp = %.10f\n", s->slowsoilp);
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
    } else if (MATCH("state", "stemp")) {
        fprintf(c->ofp, "stemp = %.10f\n", s->stemp);
        *match = TRUE;
    } else if (MATCH("state", "stempimm")) {
        fprintf(c->ofp, "stempimm = %.10f\n", s->stempimm);
        *match = TRUE;
    } else if (MATCH("state", "stempmob")) {
        fprintf(c->ofp, "stempmob = %.10f\n", s->stempmob);
        *match = TRUE;
    } else if (MATCH("state", "structsoil")) {
        fprintf(c->ofp, "structsoil = %.10f\n", s->structsoil);
        *match = TRUE;
    } else if (MATCH("state", "structsoiln")) {
        fprintf(c->ofp, "structsoiln = %.10f\n", s->structsoiln);
        *match = TRUE;
    } else if (MATCH("state", "structsoilp")) {
        fprintf(c->ofp, "structsoilp = %.10f\n", s->structsoilp);
        *match = TRUE;
    } else if (MATCH("state", "structsurf")) {
        fprintf(c->ofp, "structsurf = %.10f\n", s->structsurf);
        *match = TRUE;
    } else if (MATCH("state", "structsurfn")) {
        fprintf(c->ofp, "structsurfn = %.10f\n", s->structsurfn);
        *match = TRUE;
    } else if (MATCH("state", "structsurfp")) {
        fprintf(c->ofp, "structsurfp = %.10f\n", s->structsurfp);
        *match = TRUE;
    }

    return (1);
}
