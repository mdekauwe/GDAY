/* ============================================================================
* Soil C, N and P flows into 4 litter pools (structural and metabolic, both
* above and belowground) and 3 SOM pools (Active, slow and passive). Soil P
* flows into 5 inorganic pools (parent, lab, sorb, ssorb, and occluded).
* In essence the CENTURY model.
*
* Active pool -> soil microbes & microbial products, turnover time of mths-yrs.
* Slow pool -> resistant plant material, turnover time of 20-50 yrs.
* Passive pool -> very resistant to decomp, turnover time of > 400 yrs.
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   14.08.2016
*
* =========================================================================== */
#include "soils.h"

void calculate_csoil_flows(control *c, fluxes *f, params *p, state *s,
                           double tsoil, int doy) {
    double lnleaf, lnroot;
    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);

    /* need to store grazing flag. Allows us to switch on the annual
       grazing event, but turn it off for every other day of the year.  */
    int cntrl_grazing = c->grazing;
    if (c->grazing == 2 && p->disturbance_doy == doy+1) {
        c->grazing = TRUE;
    }

    f->tfac_soil_decomp = calc_soil_temp_factor(tsoil);

    /* calculate model decay rates */
    calculate_decay_rates(f, p, s);

    /*
     * plant litter inputs to the metabolic and structural pools determined
     * by ratio of lignin/N ratio
     */
    lnleaf = calc_ligin_nratio_leaves(c, f, p);
    lnroot = calc_ligin_nratio_fine_roots(c, f, p);
    p->fmleaf = metafract(lnleaf);
    p->fmroot = metafract(lnroot);

    /* input from faeces */
    flux_from_grazers(c, f, p);
    partition_plant_litter(f, p);
    cfluxes_from_structural_pool(f, p, s);
    cfluxes_from_metabolic_pool(f, p, s);
    cfluxes_from_active_pool(f, p, s, frac_microb_resp);
    cfluxes_from_slow_pool(f, p, s);
    cfluxes_from_passive_pool(f, p, s);
    calculate_soil_respiration(c, f, p, s);

    /* update the C pools */
    calculate_cpools(f, s);

    /* calculate NEP */
    f->nep = f->npp - f->hetero_resp - f->ceaten * (1.0 - p->fracfaeces);

    /* save fluxes for NCEAS output */
    f->co2_rel_from_surf_struct_litter = f->co2_to_air[0];
    f->co2_rel_from_soil_struct_litter = f->co2_to_air[1];
    f->co2_rel_from_surf_metab_litter = f->co2_to_air[2];
    f->co2_rel_from_soil_metab_litter = f->co2_to_air[3];
    f->co2_rel_from_active_pool = f->co2_to_air[4];
    f->co2_rel_from_slow_pool = f->co2_to_air[5];
    f->co2_rel_from_passive_pool = f->co2_to_air[6];

    /* switch off grazing if this was just activated as an annual event */
    c->grazing = cntrl_grazing;

    if (c->exudation && c->alloc_model != GRASSES) {
        calc_root_exudation_uptake_of_C(f, p, s);
    }

    return;
}

void calc_root_exudation_uptake_of_C(fluxes *f, params *p, state *s) {
    /* The amount of C which enters the active pool varies according to the
    CUE of SOM in response to root exudation (REXCUE). REXCUE determines
    the fraction of REXC that enters the active pool as C. The remaining
    flux is respired.


    REXCUE determines which fraction of REXC enters the active pool as C
    (delta_Cact). The remaining fraction of REXC is respired as CO2.
    */
    double active_CN, rex_NC, C_to_active_pool;

    active_CN = s->activesoil / s->activesoiln;

    if (p->root_exu_CUE < -0.5) {
        /*
            flexible cue
             - The constraint of 0.3<=REXCUE<=0.6 is based on observations of
               the physical limits of microbes
        */

        if (float_eq(f->root_exc, 0.0)) {
            rex_NC = 0.0;
        } else {
            rex_NC = f->root_exn / f->root_exc;
        }
        f->rexc_cue = MAX(0.3, MIN(0.6, rex_NC * active_CN));
    } else {
        f->rexc_cue = p->root_exu_CUE;
    }

    C_to_active_pool = f->root_exc * f->rexc_cue;
    s->activesoil += C_to_active_pool;

    /* Update respiration fluxes. */

    /*
    ** CUE of microbial rhizodeposition uptake is constant, so the fraction
    ** of the rhizodeposition will be used for immediate respiration
    */
    f->co2_released_exud = (1.0 - f->rexc_cue) * f->root_exc;
    f->hetero_resp += f->co2_released_exud;

    return;
}

void calculate_decay_rates(fluxes *f, params *p, state *s) {
    /* Model decay rates - decomposition rates have a strong temperature
    and moisture dependency. Note same temperature is assumed for all 3
    SOM pools, found by Knorr et al (2005) to be untrue. N mineralisation
    depends on top soil moisture (most variable) (Connell et al. 1995)

    References:
    -----------
    Knorr et al. (2005) Nature, 433, 298-301.
    Connell et al. (1995) Biol. Fert. Soils, 20, 213-220.
    */
    double adfac, soil_text, lignin_cont_leaf, lignin_cont_root;

    /* abiotic decomposition factor - impact of soil moisture
       and soil temperature on microbial activity */
    adfac = s->wtfac_root * f->tfac_soil_decomp;

    /*  Effect of soil texture (silt + clay content) on active SOM turnover
        -> higher turnover for sandy soils */
    soil_text = 1.0 - (0.75 * p->finesoil);

    /* Impact of lignin content */
    lignin_cont_leaf = exp(-3.0 * p->ligshoot);
    lignin_cont_root = exp(-3.0 * p->ligroot);

    /* decay rate of surface structural pool */
    p->decayrate[0] = p->kdec1 * lignin_cont_leaf * adfac;

    /* decay rate of surface metabolic pool */
    p->decayrate[1] = p->kdec2 * adfac;

    /* decay rate of soil structural pool */
    p->decayrate[2] = p->kdec3 * lignin_cont_root * adfac;

    /* decay rate of soil metabolic pool */
    p->decayrate[3] = p->kdec4 * adfac;

    /* decay rate of active pool */
    p->decayrate[4] = p->kdec5 * soil_text * adfac;

    /* decay rate of slow pool */
    p->decayrate[5] = p->kdec6 * adfac;

    /* decay rate of passive pool */
    p->decayrate[6] = p->kdec7 * adfac;

    return;
}

double calc_soil_temp_factor(double tsoil) {
    /*
    Soil-temperature activity factor (A9). Fit to Parton's fig 2a

    Parameters:
    -----------
    tsoil : double
        soil temperature (deg C)

    Returns:
    --------
    tfac : double
        soil temperature factor [degC]

    */
    double tfac;
    if (tsoil > 0.0) {
        tfac = MAX(0.0, 0.0326 + 0.00351 * pow(tsoil, 1.652) - \
                        pow((tsoil / 41.748), 7.19));
    } else {
        /* negative number cannot be raised to a fractional power
           number would need to be complex */
        tfac = 0.0;
    }

    return (tfac);
}

void flux_from_grazers(control *c, fluxes *f, params *p) {
    /*  Input from faeces */
    if (c->grazing) {
        p->fmfaeces = metafract(p->ligfaeces * p->faecescn / p->cfracts);
        f->faecesc = f->ceaten * p->fracfaeces;
    } else {
        p->fmfaeces = 0.0;
        f->faecesc = 0.0;
    }

    return;
}

double calc_ligin_nratio_leaves(control *c, fluxes *f, params *p) {
    /* Estimate Lignin/N ratio, as this dictates the how plant litter is
    seperated between metabolic and structural pools.

    Returns:
    --------
    lnleaf : float
        lignin:N ratio of leaf

    */
    double lnleaf, nc_leaf_litter;

    nc_leaf_litter = ratio_of_litternc_to_live_leafnc(c, f, p);

    if (float_eq(nc_leaf_litter, 0.0)) {
        /* catch divide by zero if we have no leaves */
        lnleaf = 0.0;
    } else {
        lnleaf = p->ligshoot / p->cfracts / nc_leaf_litter;
    }

    return (lnleaf);

}

double calc_ligin_nratio_fine_roots(control *c, fluxes *f, params *p) {
    /* Estimate Lignin/N ratio, as this dictates the how plant litter is
    seperated between metabolic and structural pools.

    Returns:
    --------
    lnroot : float
        lignin:N ratio of fine root
    */
    double lnroot, nc_root_litter;

    nc_root_litter = ratio_of_litternc_to_live_rootnc(c, f, p);


    if (float_eq(nc_root_litter, 0.0)) {
        /* catch divide by zero if we have no roots */
        lnroot = 0.0;
    } else {
        lnroot = p->ligroot / p->cfracts / nc_root_litter;
    }

    return (lnroot);
}

double ratio_of_litternc_to_live_leafnc(control *c, fluxes *f, params *p) {
    /* ratio of litter N:C to live leaf N:C

    Returns:
    --------
    nc_leaf_litter : float
        N:C ratio of litter to foliage

    */
    double nc_leaf_litter;

    if (c->use_eff_nc){
        nc_leaf_litter = p->liteffnc * (1.0 - p->fretrans);
    } else {
        if (float_eq(f->deadleaves, 0.0)){
            nc_leaf_litter = 0.0;
        } else {
            nc_leaf_litter = f->deadleafn / f->deadleaves;
        }
    }
    return (nc_leaf_litter);
}

double ratio_of_litternc_to_live_rootnc(control *c, fluxes *f, params *p) {
    /* ratio of litter N:C to live root N:C

    Returns:
    --------
    nc_root_litter : float
        N:C ratio of litter to live root

    */
    double nc_root_litter;

    if (c->use_eff_nc){
        nc_root_litter = p->liteffnc * p->ncrfac *  (1.0 - p->rretrans);
    } else {
        if (float_eq(f->deadroots, 0.0)){
            nc_root_litter = 0.0;
        } else{
            nc_root_litter = f->deadrootn / f->deadroots;
        }
    }

    return (nc_root_litter);
}

double metafract(double lig2n) {
    /* Calculate what fraction of the litter will be partitioned to the
    metabolic pool which is given by the lignin:N ratio.

    Parameters:
    -----------
    lig2n : float
        lignin to N ratio

    Returns:
    --------
    metabolic fraction : float
        partitioned fraction to metabolic pool [must be positive]
    */

    /* Original implementation based on Parton et al. */
    return (MAX(0.0, 0.85 - (0.018 * lig2n)));
}


void partition_plant_litter(fluxes *f, params *p) {
    /* Partition litter from the plant (surface) and roots into metabolic
    and structural pools  */

    double leaf_material, wood_material, faeces_material;
    /*
     * Surface (leaves, branches, stem) Litter
     */

    /* ...to the structural pool*/
    leaf_material = f->deadleaves * (1.0 - p->fmleaf);
    wood_material = f->deadbranch + f->deadstems;
    faeces_material = f->faecesc * (1.0 - p->fmfaeces);
    f->surf_struct_litter = leaf_material + wood_material + faeces_material;

    /* ...to the metabolic pool */
    f->surf_metab_litter = f->deadleaves * p->fmleaf + f->faecesc * p->fmfaeces;

    /*
    ** Root Litter
    */

    /* ...to the structural pool */
    f->soil_struct_litter = f->deadroots * (1.0 - p->fmroot) + f->deadcroots;

    /* ...to the metabolic pool */
    f->soil_metab_litter = f->deadroots * p->fmroot;

    return;
}

void cfluxes_from_structural_pool(fluxes *f, params *p, state *s) {

    /* Send structural c fluxes to other SOM pools */

    double structout_surf = s->structsurf * p->decayrate[0];
    double structout_soil = s->structsoil * p->decayrate[2];

    /* C flux surface structural pool -> slow pool */
    f->surf_struct_to_slow = structout_surf * p->ligshoot * 0.7;

    /* C flux surface structural pool -> active pool */
    f->surf_struct_to_active = structout_surf * (1.0 - p->ligshoot) * 0.55;

    /* C flux soil structural pool -> slow pool */
    f->soil_struct_to_slow = structout_soil * p->ligroot * 0.7;

    /* soil structural pool -> active pool */
    f->soil_struct_to_active = structout_soil * (1.0 - p->ligroot) * 0.45;

    /* Respiration fluxes */

    /* CO2 lost during transfer of structural C to the slow pool */
    f->co2_to_air[0] = (structout_surf *
                        (p->ligshoot * 0.3 + (1.0 - p->ligshoot) * 0.45));

    /* CO2 lost during transfer structural C  to the active pool */
    f->co2_to_air[1] = (structout_soil *
                        (p->ligroot * 0.3 + (1.0 - p->ligroot) * 0.55));

    return;
}

void cfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {
    /* Send C from metabolic pools to other SOM pools */

    /* C flux surface metabolic pool -> active pool */
    f->surf_metab_to_active = s->metabsurf * p->decayrate[1] * 0.45;

    /* C flux soil metabolic pool  -> active pool */
    f->soil_metab_to_active = s->metabsoil * p->decayrate[3] * 0.45;

    /* Respiration fluxes */
    f->co2_to_air[2] = s->metabsurf * p->decayrate[1] * 0.55;
    f->co2_to_air[3] = s->metabsoil * p->decayrate[3] * 0.55;

    return;
}

void cfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {
    /* Send C fluxes from active pool to other SOM pools */

    double activeout = s->activesoil * p->decayrate[4];

    /* C flux active pool -> slow pool */
    f->active_to_slow = activeout * (1.0 - frac_microb_resp - 0.004);

    /* (Parton 1993)
    f->active_to_slow = (activeout * (1.0 - frac_microb_resp - 0.003 -
                         0.032 * Claysoil));
    */

    /* C flux active pool -> passive pool */
    f->active_to_passive = activeout * 0.004;

    /* Respiration fluxes */
    f->co2_to_air[4] = activeout * frac_microb_resp;

    return;
}

void cfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
    /* Send C fluxes from slow pool to other SOM pools */

    double slowout = s->slowsoil * p->decayrate[5];

    /* C flux slow pool -> active pool */
    f->slow_to_active = slowout * 0.42;

    /* slow pool -> passive pool */
    f->slow_to_passive = slowout * 0.03;

    /* Respiration fluxes */
    f->co2_to_air[5] = slowout * 0.55;

    return;
}

void cfluxes_from_passive_pool(fluxes *f, params *p, state *s) {

    /* C flux passive pool -> active pool */
    f->passive_to_active = s->passivesoil * p->decayrate[6] * 0.45;

    /* Respiration fluxes */
    f->co2_to_air[6] = s->passivesoil * p->decayrate[6] * 0.55;

    return;
}

void calculate_soil_respiration(control *c, fluxes *f, params *p, state *s) {
    /* calculate the total soil respiration (heterotrophic) flux, i.e.
    the amount of CO2 released back to the atmosphere */

    /* total CO2 production */
    f->hetero_resp = (f->co2_to_air[0] + f->co2_to_air[1] + f->co2_to_air[2] +
                      f->co2_to_air[3] + f->co2_to_air[4] + f->co2_to_air[5] +
                      f->co2_to_air[6]);

    /* insert following line so value of respiration obeys c conservation if
       assuming a fixed passive pool */
    if (c->passiveconst == TRUE) {
        f->hetero_resp = (f->hetero_resp + f->active_to_passive +
                          f->slow_to_passive - s->passivesoil *
                          p->decayrate[6]);
    }

    return;
}

void calculate_cpools(fluxes *f, state *s) {
    /* Calculate new soil carbon pools. */

    /* Update pools */
    s->structsurf += (f->surf_struct_litter -
                     (f->surf_struct_to_slow + f->surf_struct_to_active +
                      f->co2_to_air[0]));

    s->structsoil += (f->soil_struct_litter -
                     (f->soil_struct_to_slow + f->soil_struct_to_active +
                      f->co2_to_air[1]));

    s->metabsurf += (f->surf_metab_litter -
                     (f->surf_metab_to_active + f->co2_to_air[2]));

    s->metabsoil += (f->soil_metab_litter -
                     (f->soil_metab_to_active + f->co2_to_air[3]));

    /* store the C SOM fluxes for Nitrogen/Phosphorus calculations */
    f->c_into_active = (f->surf_struct_to_active + f->soil_struct_to_active +
                        f->surf_metab_to_active + f->soil_metab_to_active +
                        f->slow_to_active + f->passive_to_active);

    f->c_into_slow = (f->surf_struct_to_slow + f->soil_struct_to_slow +
                      f->active_to_slow);

    f->c_into_passive = f->active_to_passive + f->slow_to_passive;

    s->activesoil += (f->c_into_active -
                      (f->active_to_slow + f->active_to_passive +
                       f->co2_to_air[4]));

    s->slowsoil += (f->c_into_slow -
                    (f->slow_to_active + f->slow_to_passive +
                     f->co2_to_air[5]));

    s->passivesoil += (f->c_into_passive -
                        (f->passive_to_active + f->co2_to_air[6]));

    /*
      When nothing is being added to the metabolic pools, there is the
      potential scenario with the way the model works for tiny bits to be
      removed with each timestep. Effectively with time this value which is
      zero can end up becoming zero but to a silly decimal place
    */
    precision_control_soil_c(f, s);

    return;
}

void precision_control_soil_c(fluxes *f, state *s) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-08, excess;

    /* C & N state variables */
    if (s->metabsurf < tolerance) {
        excess = s->metabsurf;
        f->surf_metab_to_active = excess * 0.45;
        f->co2_to_air[2] = excess * 0.55;
        s->metabsurf = 0.0;
    }

    if (s->metabsoil < tolerance) {
        excess = s->metabsoil;
        f->soil_metab_to_active = excess * 0.45;
        f->co2_to_air[3] = excess * 0.55;
        s->metabsoil = 0.0;
    }

    return;
}


void calculate_nsoil_flows(control *c, fluxes *f, params *p, state *s,
                           int doy) {

    /* need to store grazing flag. Allows us to switch on the annual
       grazing event, but turn it off for every other day of the year.  */
    int cntrl_grazing = c->grazing;
    if (c->grazing == 2 && p->disturbance_doy == doy+1) {
        c->grazing = TRUE;
    }

    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double nsurf, nsoil, active_nc_slope, slow_nc_slope, passive_nc_slope;

    grazer_inputs(c, f, p);
    n_inputs_from_plant_litter(f, p, &nsurf, &nsoil);
    partition_plant_litter_n(c, f, p, nsurf, nsoil);

    /* SOM nitrogen effluxes.  These are assumed to have the source n:c
       ratio prior to the increase of N:C due to co2 evolution. */
    nfluxes_from_structural_pools(f, p, s);
    nfluxes_from_metabolic_pool(f, p, s);
    nfluxes_from_active_pool(f, p, s, frac_microb_resp);
    nfluxes_from_slow_pool(f, p, s);
    nfluxes_from_passive_pool(f, p, s);

    /* gross N mineralisation */
    calculate_n_mineralisation(f);

    /* calculate N immobilisation */
    calculate_n_immobilisation(f, p, s, &(f->nimmob), &active_nc_slope,
                               &slow_nc_slope, &passive_nc_slope);

    /* calculate N net mineralisation */
    calc_n_net_mineralisation(f);

    if (c->exudation && c->alloc_model != GRASSES) {
        calc_root_exudation_uptake_of_N(f, s);
    }

    if (c->adjust_rtslow) {
        adjust_residence_time_of_slow_pool(f, p);
    } else {
        /* Need to correct units of rate constant */
        f->rtslow = 1.0 / (p->kdec6 * NDAYS_IN_YR);
    }

    /* Update model soil N pools */
    calculate_npools(c, f, p, s, active_nc_slope, slow_nc_slope,
                     passive_nc_slope);

    /* switch off grazing if this was just activated as an annual event */
    c->grazing = cntrl_grazing;

    return;
}

void calc_root_exudation_uptake_of_N(fluxes *f, state *s) {
    /* When N mineralisation is large enough to allow a small amount of N
    immobilisation, the amount of N which enters the active pool is
    calculated according to REXC divided by the CN of the active pool. When
    exudation enters the active pool, the CN ratio of the exudates drops
    from REXC/REXN to the CN of the active pool. Which is consistent with
    the CENTURY framework, where C flows between pools lead to either
    mineralisation (N gain) or immobilisation (N loss) due to differences
    in the CN ratio of the outgoing and incoming pools.

    The amount of N added to the active pool is independent of the CUE of
    the microbial pool in response to root exudation (REXCUE).
    */
    double N_available, active_NC, delta_Nact, N_miss, N_to_active_pool;

    N_available = s->inorgn + (f->ninflow + f->nurine + f->nmineralisation -
                               f->nloss - f->nuptake);

    active_NC = s->activesoiln / s->activesoil;
    delta_Nact = f->root_exc * f->rexc_cue * active_NC;

    /*
    ** Demand for N from exudation to meet the C:N ratio of the active pool,
    ** given the amount of N you add.
    */
    N_miss = delta_Nact - f->root_exn;

    if (N_miss <= 0.0) {
        /*
        ** Root exudation includes more N than is needed by the microbes, the
        ** excess is mineralised
        */
        f->nmineralisation -= N_miss;
        N_to_active_pool = f->root_exn + N_miss;
    } else {
        /*
        ** Not enough N in the soil to meet demand, so we are providing all
        ** the N we have, which means that the C:N ratio of the active pool
        ** changes.
        */
        if (N_miss > N_available) {
            N_to_active_pool = f->root_exn + N_available;
            f->nmineralisation -= N_available;
        } else {
            /*
            ** Enough N to meet demand, so takes N from the mineralisation
            ** and the active pool maintains the same C:N ratio.
            */
            N_to_active_pool = f->root_exn + N_miss;
            f->nmineralisation -= N_miss;
        }
    }

    /* update active pool */
    s->activesoiln += N_to_active_pool;

    return;
}

void adjust_residence_time_of_slow_pool(fluxes *f, params *p) {
    /* Priming simulations the residence time of the slow pool is flexible,
    as the flux out of the active pool (factive) increases the residence
    time of the slow pool decreases.
    */
    double rt_slow_pool;

    /* total flux out of the factive pool */
    f->factive = (f->active_to_slow + f->active_to_passive + \
                  f->co2_to_air[4] + f->co2_released_exud);

    if (float_eq(f->factive, 0.0)) {
        /* Need to correct units of rate constant */
        rt_slow_pool = 1.0 / (p->kdec6 * NDAYS_IN_YR);
    } else {
        rt_slow_pool = (1.0 / p->prime_y) / \
                        MAX(0.3, (f->factive / (f->factive + p->prime_z)));

        /* GDAY uses decay rates rather than residence times... */
        p->kdec6 = 1.0 / rt_slow_pool;

        /* rate constant needs to be per day inside GDAY */
        p->kdec6 /= NDAYS_IN_YR;

    }

    /* Save for outputting purposes only */
    f->rtslow = rt_slow_pool;

    return;
}

void grazer_inputs(control *c, fluxes *f, params *p) {
    /* Grazer inputs from faeces and urine, flux detd by faeces c:n */
    double arg;

    if (c->grazing)
        p->faecesn = f->faecesc / p->faecescn;
    else
        p->faecesn = 0.0;

    /* make sure faecesn <= total n input to soil from grazing */
    arg = f->neaten * p->fractosoil;
    if (p->faecesn > arg)
        p->faecesn = f->neaten * p->fractosoil;

    /* urine=total-faeces */
    if (c->grazing)
        f->nurine = f->neaten * p->fractosoil - p->faecesn;
    else
        f->nurine = 0.0;

    if (f->nurine < 0.0)
        f->nurine = 0.0;

    return;
}

void n_inputs_from_plant_litter(fluxes *f, params *p, double *nsurf,
                              double *nsoil) {
    /* inputs from plant litter.

    surface and soil pools are independent. Structural input flux n:c can
    be either constant or a fixed fraction of metabolic input flux.

    Returns:
    --------
    nsurf : float
        N input from surface pool
    nsoil : float
        N input from soil pool
    */

    /* surface and soil inputs (faeces n goes to abovgrd litter pools) */
    *nsurf = f->deadleafn + f->deadbranchn + f->deadstemn + p->faecesn;
    *nsoil = f->deadrootn + f->deadcrootn;

    return;
}

void partition_plant_litter_n(control *c, fluxes *f, params *p, double nsurf,
                              double nsoil) {
    /* Partition litter N from the plant (surface) and roots into metabolic
    and structural pools

    Parameters:
    -----------
    nsurf : float
        N input from surface pool
    nsoil : float
        N input from soil pool
    */

    double c_surf_struct_litter, c_soil_struct_litter;

    /* constant structural input n:c as per century */
    if (c->strfloat == 1) {

        /* structural input n:c is a fraction of metabolic */
        c_surf_struct_litter = (f->surf_struct_litter * p->structrat +
                                f->surf_metab_litter);

        if (float_eq(c_surf_struct_litter, 0.0))
             f->n_surf_struct_litter = 0.0;
        else
             f->n_surf_struct_litter = (nsurf * f->surf_struct_litter *
                                        p->structrat / c_surf_struct_litter);

        c_soil_struct_litter = (f->soil_struct_litter * p->structrat +
                                f->soil_metab_litter);

        if (float_eq(c_soil_struct_litter, 0.0))
            f->n_soil_struct_litter = 0.0;
        else
            f->n_soil_struct_litter = (nsurf * f->soil_struct_litter *
                                       p->structrat / c_soil_struct_litter);
    } else {

        /* n flux -> surface structural pool */
        f->n_surf_struct_litter = f->surf_struct_litter / p->structcn;

        /* n flux -> soil structural pool */
        f->n_soil_struct_litter = f->soil_struct_litter / p->structcn;

        /* if not enough N for structural, all available N goes to structural */
        if (f->n_surf_struct_litter > nsurf)
             f->n_surf_struct_litter = nsurf;
        if (f->n_soil_struct_litter > nsoil)
            f->n_soil_struct_litter = nsoil;
    }


    /* remaining N goes to metabolic pools */
    f->n_surf_metab_litter = nsurf - f->n_surf_struct_litter;
    f->n_soil_metab_litter = nsoil - f->n_soil_struct_litter;

    return;
}

void nfluxes_from_structural_pools(fluxes *f, params *p, state *s) {
    /* from structural pool */
    double sigwt;
    double structout_surf = s->structsurfn * p->decayrate[0];
    double structout_soil = s->structsoiln * p->decayrate[2];

    sigwt = structout_surf / (p->ligshoot * 0.7 + (1.0 - p->ligshoot) * 0.55);

    /* N flux from surface structural pool -> slow pool */
    f->n_surf_struct_to_slow = sigwt * p->ligshoot * 0.7;

    /* N flux surface structural pool -> active pool */
    f->n_surf_struct_to_active = sigwt * (1.0 - p->ligshoot) * 0.55;

    sigwt = structout_soil / (p->ligroot * 0.7 + (1. - p->ligroot) * 0.45);


    /* N flux from soil structural pool -> slow pool */
    f->n_soil_struct_to_slow = sigwt * p->ligroot * 0.7;

    /* N flux from soil structural pool -> active pool */
    f->n_soil_struct_to_active = sigwt * (1.0 - p->ligroot) * 0.45;

    return;
}

void nfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {

    /* N flux surface metabolic pool -> active pool */
    f->n_surf_metab_to_active = s->metabsurfn * p->decayrate[1];

    /* N flux soil metabolic pool  -> active pool */
    f->n_soil_metab_to_active = s->metabsoiln * p->decayrate[3];

    return;
}

void nfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {

    double activeout, sigwt;
    /* N fluxes from active pool */
    activeout = s->activesoiln * p->decayrate[4];
    sigwt = activeout / (1.0 - frac_microb_resp);

    /* N flux active pool -> slow pool */
    f->n_active_to_slow = sigwt * (1.0 - frac_microb_resp - 0.004);

    /* N flux active pool -> passive pool */
    f->n_active_to_passive = sigwt * 0.004;

    return;
}

void nfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
    /* N fluxes from slow pools */

    double slowout = s->slowsoiln * p->decayrate[5];
    double sigwt = slowout / 0.45;

    /* C flux slow pool -> active pool */
    f->n_slow_to_active = sigwt * 0.42;

    /* slow pool -> passive pool */
    f->n_slow_to_passive = sigwt * 0.03;

    return;
}

void nfluxes_from_passive_pool(fluxes *f, params *p, state *s) {
    /* N fluxes from passive pool */

    /* C flux passive pool -> active pool */
    f->n_passive_to_active = s->passivesoiln * p->decayrate[6];

    return;
}

void calculate_n_mineralisation(fluxes *f) {
    /* N gross mineralisation rate is given by the excess of N outflows
    over inflows. Nitrogen mineralisation is the process by which organic
    N is converted to plant available inorganic N, i.e. microbes decompose
    organic N from organic matter to ammonia (NH3) and ammonium (NH4),
    called ammonification.

    Returns:
    --------
    value : float
        Gross N mineralisation
    */
    f->ngross =  (f->n_surf_struct_to_slow + f->n_surf_struct_to_active +
                  f->n_soil_struct_to_slow + f->n_soil_struct_to_active +
                  f->n_surf_metab_to_active + f->n_soil_metab_to_active +
                  f->n_active_to_slow + f->n_active_to_passive +
                  f->n_slow_to_active + f->n_slow_to_passive +
                  f->n_passive_to_active);
    return;
}

void calculate_n_immobilisation(fluxes *f, params *p, state *s, double *nimmob,
                                double *active_nc_slope, double *slow_nc_slope,
                                double *passive_nc_slope) {
    /* N immobilised in new soil organic matter, the reverse of
    mineralisation. Micro-organisms in the soil compete with plants for N.
    Immobilisation is the process by which nitrate and ammonium are taken up
    by the soil organisms and thus become unavailable to the plant
    (->organic N).

    When C:N ratio is high the microorganisms need more nitrogen from
    the soil to decompose the carbon in organic materials. This N will be
    immobilised until these microorganisms die and the nitrogen is
    released.

    General equation for new soil N:C ratio vs Nmin, expressed as linear
    equation passing through point Nmin0, actncmin (etc). Values can be
    Nmin0=0, Actnc0=Actncmin

    if Nmin < Nmincrit:
        New soil N:C = soil N:C (when Nmin=0) + slope * Nmin

    if Nmin > Nmincrit
        New soil N:C = max soil N:C

    NB N:C ratio of new passive SOM can change even if assume Passiveconst

    Returns:
    --------
    nimob : float
        N immobilsed
    */
    double nmin, arg1, arg2, arg3, numer1, numer2, denom;

    /* N:C new SOM - active, slow and passive */
    *active_nc_slope = calculate_nc_slope(p, p->actncmax, p->actncmin);
    *slow_nc_slope = calculate_nc_slope(p, p->slowncmax, p->slowncmin);
    *passive_nc_slope = calculate_nc_slope(p, p->passncmax, p->passncmin);

    /* convert units */
    nmin = p->nmin0 / M2_AS_HA * G_AS_TONNES;

    arg1 = (p->passncmin - *passive_nc_slope * nmin) * f->c_into_passive;
    arg2 = (p->slowncmin - *slow_nc_slope * nmin) * f->c_into_slow;
    arg3 = f->c_into_active * (p->actncmin - *active_nc_slope * nmin);
    numer1 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * p->passncmax;
    arg2 = f->c_into_slow * p->slowncmax;
    arg3 = f->c_into_active * p->actncmax;
    numer2 = arg1 + arg2 + arg3;

    arg1 = f->c_into_passive * *passive_nc_slope;
    arg2 = f->c_into_slow * *slow_nc_slope;
    arg3 = f->c_into_active * *active_nc_slope;
    denom = arg1 + arg2 + arg3;

    /* evaluate N immobilisation in new SOM */
    *nimmob = numer1 + denom * s->inorgn;
    if (*nimmob > numer2)
        *nimmob = numer2;

    return;
}


void calc_n_net_mineralisation(fluxes *f) {
    /* N Net mineralisation from microbial activity */
    f->nmineralisation = f->ngross - f->nimmob + f->nlittrelease;

    return;
}

double calculate_nc_slope(params *p, double ncmax, double ncmin) {
    /* Returns N:C ratio of the mineral pool slope

    based on fig 4 of Parton et al 1993. Standard slow pool C:N is different
    to the values in Parton. Bill quotes 12-20, whereas McMurtrie et al '01
    use 10-40.

    Parameters
    ----------
    ncmax : float
        SOM pools maximum N:C
    ncmin: float
        SOM pools minimum N:C

    Returns:
    --------
    value : float
        SOM pool N:C ratio
    */
    double arg1, arg2, conv;

    arg1 = ncmax - ncmin;
    arg2 = p->nmincrit - p->nmin0;
    conv = M2_AS_HA / G_AS_TONNES;

    return (arg1 / arg2 * conv);
}

void calculate_npools(control *c, fluxes *f, params *p, state *s,
                      double active_nc_slope, double slow_nc_slope,
                      double passive_nc_slope) {
    /*
    Update N pools in the soil

    Parameters
    ----------
    active_nc_slope : float
        active NC slope
    slow_nc_slope: float
        slow NC slope
    passive_nc_slope : float
        passive NC slope

    */
    double n_into_active, n_out_of_active, n_into_slow, n_out_of_slow,
           n_into_passive, n_out_of_passive, arg, active_nc, fixn, slow_nc,
           pass_nc;

    /*
        net N release implied by separation of litter into structural
        & metabolic. The following pools only fix or release N at their
        limiting n:c values.
    */

    /* N released or fixed from the N inorganic pool is incremented with
       each call to nc_limit and stored in f->nlittrelease */
    f->nlittrelease = 0.0;

    s->structsurfn += (f->n_surf_struct_litter -
                        (f->n_surf_struct_to_slow +
                         f->n_surf_struct_to_active));

    s->structsoiln += (f->n_soil_struct_litter -
                       (f->n_soil_struct_to_slow + f->n_soil_struct_to_active));

    if (c->strfloat == 0) {
        s->structsurfn += nc_limit(f, s->structsurf, s->structsurfn,
                                   1.0/p->structcn, 1.0/p->structcn);
        s->structsoiln += nc_limit(f, s->structsoil, s->structsoiln,
                                   1.0/p->structcn, 1.0/p->structcn);
    }

    s->metabsurfn += f->n_surf_metab_litter - f->n_surf_metab_to_active;
    s->metabsurfn += nc_limit(f, s->metabsurf, s->metabsurfn,1.0/25.0, 1.0/10.0);
    s->metabsoiln += (f->n_soil_metab_litter - f->n_soil_metab_to_active);
    s->metabsoiln += nc_limit(f, s->metabsoil, s->metabsoiln, 1.0/25.0,
                              1.0/10.0);

    /* When nothing is being added to the metabolic pools, there is the
       potential scenario with the way the model works for tiny bits to be
       removed with each timestep. Effectively with time this value which is
       zero can end up becoming zero but to a silly decimal place */
    precision_control_soil_n(f, s);

    /* Update SOM pools */
    n_into_active = (f->n_surf_struct_to_active + f->n_soil_struct_to_active +
                     f->n_surf_metab_to_active + f->n_soil_metab_to_active +
                     f->n_slow_to_active + f->n_passive_to_active);

    n_out_of_active = f->n_active_to_slow + f->n_active_to_passive;

    n_into_slow = (f->n_surf_struct_to_slow + f->n_soil_struct_to_slow +
                   f->n_active_to_slow);

    n_out_of_slow = f->n_slow_to_active + f->n_slow_to_passive;
    n_into_passive = f->n_active_to_passive + f->n_slow_to_passive;
    n_out_of_passive = f->n_passive_to_active;

    /* N:C of the SOM pools increases linearly btw prescribed min and max
       values as the Nconc of the soil increases. */
    arg = s->inorgn - p->nmin0 / M2_AS_HA * G_AS_TONNES;

    /* active */
    active_nc = p->actncmin + active_nc_slope * arg;
    if (active_nc > p->actncmax)
        active_nc = p->actncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_active, n_into_active, active_nc);
    s->activesoiln += n_into_active + fixn - n_out_of_active;

    /* slow */
    slow_nc = p->slowncmin + slow_nc_slope * arg;
    if (slow_nc > p->slowncmax)
        slow_nc = p->slowncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_slow, n_into_slow, slow_nc);
    s->slowsoiln += n_into_slow + fixn - n_out_of_slow;

    /* passive, update passive pool only if passiveconst=0 */
    pass_nc = p->passncmin + passive_nc_slope * arg;
    if (pass_nc > p->passncmax)
        pass_nc = p->passncmax;

    /* release N to Inorganic pool or fix N from the Inorganic pool in order
       to normalise the N:C ratio of a net flux */
    fixn = nc_flux(f->c_into_passive, n_into_passive, pass_nc);
    s->passivesoiln += n_into_passive + fixn - n_out_of_passive;

    /* Daily increment of soil inorganic N pool, diff btw in and effluxes
       (grazer urine n goes directly into inorganic pool) nb inorgn may be
       unstable if rateuptake is large */
    s->inorgn += (f->ninflow + f->nurine + f->nmineralisation -
                  f->nloss - f->nuptake);

    //fprintf(stderr, "inorgn = %f\n", s->inorgn);

    /*f->nmineralisation = f->ngross - f->nimmob + f->nlittrelease;*/
    return;
}


double nc_limit(fluxes *f, double cpool, double npool, double ncmin,
                double ncmax) {
    /* Release N to 'Inorgn' pool or fix N from 'Inorgn', in order to keep
    the  N:C ratio of a litter pool within the range 'ncmin' to 'ncmax'.

    Parameters:
    -----------
    cpool : float
        various C pool (state)
    npool : float
        various N pool (state)
    ncmin : float
        maximum N:C ratio
    ncmax : float
        minimum N:C ratio

    Returns:
    --------
    fix/rel : float
        amount of N to be added/released from the inorganic pool

    */
    double rel, fix;
    double nmax = cpool * ncmax;
    double nmin = cpool * ncmin;

    if (npool > nmax) {
        /* release */
        rel = npool - nmax;
        f->nlittrelease += rel;
        return (-rel);
    } else if (npool < nmin) {
        /* fix */
        fix = nmin - npool;
        f->nlittrelease -= fix;
        return (fix);
    } else {
        return (0.0);
    }
}

double nc_flux(double cflux, double nflux, double nc_ratio) {
    /*
    Release N to Inorganic pool or fix N from the Inorganic pool in order
    to normalise the N:C ratio of a net flux

    Parameters:
    -----------
    cflux : float
        C flux into SOM pool
    nflux : float
        N flux into SOM pool
    nc_ratio : float
        preferred N:C ratio

    Returns:
        fix : float
        Returns the amount of N required to be fixed
    */

    return (cflux * nc_ratio) - nflux;
}


void precision_control_soil_n(fluxes *f, state *s) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-08, excess;

    if (s->metabsurfn < tolerance) {
        excess = s->metabsurfn;
        f->n_surf_metab_to_active = excess;
        s->metabsurfn = 0.0;
    }

    if (s->metabsoiln < tolerance) {
        excess = s->metabsoiln;
        f->n_soil_metab_to_active = excess;
        s->metabsoiln = 0.0;
    }

    return;
}

void soil_soprtion_parameters(char *soil_order, params *p) {
    //
    // Parameterize Smax and Ks parameters based on soil order;
    // Ref. Yang et al., 2016, GRL, Table S2
    //

    if (strcmp(soil_order, "andisol") == 0) {
        p->smax = 10;
    } else if (strcmp(soil_order, "gelisol") == 0) {
        p->smax = 5;
    } else if (strcmp(soil_order, "histosol") == 0) {
        p->smax = 5;
    } else if (strcmp(soil_order, "entisol") == 0) {
        p->smax = 5;
    } else if (strcmp(soil_order, "inceptisol") == 0) {
        p->smax = 5;
    } else if (strcmp(soil_order, "aridsol") == 0) {
        p->smax = 7;
    } else if (strcmp(soil_order, "vertisol") == 0) {
        p->smax = 7;
    } else if (strcmp(soil_order, "mollisol") == 0) {
        p->smax = 7;
    } else if (strcmp(soil_order, "alfisol") == 0) {
        p->smax = 7;
    } else if (strcmp(soil_order, "spodosol") == 0) {
        p->smax = 9;
    } else if (strcmp(soil_order, "ultisol") == 0) {
        p->smax = 9;
    } else if (strcmp(soil_order, "oxisol") == 0) {
        p->smax = 9;
    } else {
        fprintf(stderr, "Could not understand soil order", __LINE__);
        exit(EXIT_FAILURE);
    }

    return;
}


void calculate_psoil_flows(control *c, fluxes *f, params *p, state *s,
                           int doy) {
    /*
        need to store grazing flag. Allows us to switch on the annual
        grazing event, but turn it off for every other day of the year.
    */
    int cntrl_grazing = c->grazing;
    if (c->grazing == 2 && p->disturbance_doy == doy+1) {
        c->grazing = TRUE;
    }

    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double psurf, psoil, active_pc_slope, slow_pc_slope, passive_pc_slope;

    grazer_inputs_p(c, f, p);
    p_inputs_from_plant_litter(f, p, &psurf, &psoil);
    partition_plant_litter_p(c, f, p, psurf, psoil);

    /* SOM phosphorus effluxes.  These are assumed to have the source p:c
    ratio prior to the increase of P:C due to co2 evolution. */
    pfluxes_from_structural_pools(f, p, s);
    pfluxes_from_metabolic_pool(f, p, s);
    pfluxes_from_active_pool(f, p, s, frac_microb_resp);
    pfluxes_from_slow_pool(f, p, s);
    pfluxes_from_passive_pool(f, p, s);

    /* calculate P parent influxe to mineral P */
    calculate_p_parent_fluxes(f, p, s);

    /* gross P mineralisation */
    calculate_p_mineralisation(f);

    /* calculate P immobilisation */
    calculate_p_immobilisation(f, p, s, &(f->pimmob), &active_pc_slope,
                             &slow_pc_slope, &passive_pc_slope);

    /* calculate P net mineralisation, excluding biochemical mineralisation */
    calc_p_net_mineralisation(f);

    /* calculate P biochemical mineralisation */
    calculate_p_biochemical_mineralisation(f, p, s);

    /* SIM phosphorus dynamics */
    calculate_p_ssorb_to_sorb(s, f, p, c);
    calculate_p_ssorb_to_occ(s, f, p);
    calculate_p_sorb_to_ssorb(s, f, p);

    /* calculate P lab and sorb fluxes from gross P flux */
    calculate_p_min_fluxes(f, p, s);

    if (c->exudation && c->alloc_model != GRASSES) {
        calc_root_exudation_uptake_of_P(f, s);
    }

    if (c->adjust_rtslow) {
        adjust_residence_time_of_slow_pool(f, p);
    } else {
        /* Need to correct units of rate constant */
        f->rtslow = 1.0 / (p->kdec6 * NDAYS_IN_YR);
    }

    /* Update model soil P pools */
    calculate_ppools(c, f, p, s, active_pc_slope, slow_pc_slope,
                     passive_pc_slope);

    /* switch off grazing if this was just activated as an annual event */
    c->grazing = cntrl_grazing;

    return;
}

void calc_root_exudation_uptake_of_P(fluxes *f, state *s) {
  /* When P mineralisation is large enough to allow a small amount of P
  immobilisation, the amount of P which enters the active pool is
  calculated according to REXC divided by the CP of the active pool. When
  exudation enters the active pool, the CP ratio of the exudates drops
  from REXC/REXP to the CP of the active pool. Which is consistent with
  the CENTURY framework, where C flows between pools lead to either
  mineralisation (P gain) or immobilisation (P loss) due to differences
  in the CP ratio of the outgoing and incoming pools.

  The amount of P added to the active pool is independent of the CUE of
  the microbial pool in response to root exudation (REXCUE).
  */
  double P_available, active_PC, delta_Pact, P_miss, P_to_active_pool;

  P_available = s->inorglabp -f->ploss - f->puptake;

  active_PC = s->activesoilp / s->activesoil;
  delta_Pact = f->root_exc * f->rexc_cue * active_PC;

  /*
  ** Demand for P from exudation to meet the C:P ratio of the active pool,
  ** given the amount of P you add.
  */
  P_miss = delta_Pact - f->root_exp;

  if (P_miss <= 0.0) {
    /*
    ** Root exudation includes more P than is needed by the microbes, the
    ** excess is mineralised
    */
    f->pmineralisation -= P_miss;
    P_to_active_pool = f->root_exp + P_miss;
  } else {
    /*
    ** Not enough P in the soil to meet demand, so we are providing all
    ** the P we have, which means that the C:P ratio of the active pool
    ** changes.
    */
    if (P_miss > P_available) {
      P_to_active_pool = f->root_exp + P_available;
      f->pmineralisation -= P_available;
    } else {
      /*
      ** Enough P to meet demand, so takes P from the mineralisation
      ** and the active pool maintains the same C:P ratio.
      */
      P_to_active_pool = f->root_exp + P_miss;
      f->pmineralisation -= P_miss;
    }
  }

  /* update active pool */
  s->activesoilp += P_to_active_pool;

  return;
}

void grazer_inputs_p(control *c, fluxes *f, params *p) {
  /* Grazer inputs from faeces and urine, flux detd by faeces c:p */
  double arg;

  if (c->grazing)
    p->faecesp = f->faecesc / p->faecescp;
  else
    p->faecesp = 0.0;

  /* make sure faecesp <= total p input to soil from grazing */
  arg = f->peaten * p->fractosoilp;
  if (p->faecesp > arg)
    p->faecesp = f->peaten * p->fractosoilp;

  /* urine=total-faeces */
  if (c->grazing)
    f->purine = f->peaten * p->fractosoilp - p->faecesp;
  else
    f->purine = 0.0;

  if (f->purine < 0.0)
    f->purine = 0.0;

  return;
}

void p_inputs_from_plant_litter(fluxes *f, params *p, double *psurf,
                              double *psoil) {
  /* inputs from plant litter.

  surface and soil pools are independent. Structural input flux p:c can
  be either constant or a fixed fraction of metabolic input flux.

  Returns:
  --------
  psurf : float
  P input from surface pool
  psoil : float
  P input from soil pool
  */

  /* surface and soil inputs (faeces p goes to abovgrd litter pools) */
  *psurf = f->deadleafp + f->deadbranchp + f->deadstemp + p->faecesp;
  *psoil = f->deadrootp + f->deadcrootp;

  return;
}

void partition_plant_litter_p(control *c, fluxes *f, params *p, double psurf,
                              double psoil) {
  /* Partition litter P from the plant (surface) and roots into metabolic
  and structural pools

  Parameters:
  -----------
  psurf : float
  P input from surface pool
  psoil : float
  P input from soil pool
  */

  double c_surf_struct_litter, c_soil_struct_litter;

  /* constant structural input p:c as per century */
  if (c->strpfloat == 1) {
    /* structural input p:c is a fraction of metabolic */
    c_surf_struct_litter = (f->surf_struct_litter * p->structratp +
    f->surf_metab_litter);

    if (float_eq(c_surf_struct_litter, 0.0))
      f->p_surf_struct_litter = 0.0;
    else
      f->p_surf_struct_litter = (psurf * f->surf_struct_litter *
        p->structratp / c_surf_struct_litter);

      c_soil_struct_litter = (f->soil_struct_litter * p->structratp +
        f->soil_metab_litter);

    if (float_eq(c_soil_struct_litter, 0.0))
      f->p_soil_struct_litter = 0.0;
    else
      f->p_soil_struct_litter = (psurf * f->soil_struct_litter *
        p->structratp / c_soil_struct_litter);
  } else {

    /* p flux -> surface structural pool */
    f->p_surf_struct_litter = f->surf_struct_litter / p->structcp;

    /* p flux -> soil structural pool */
    f->p_soil_struct_litter = f->soil_struct_litter / p->structcp;

    /* if not enough P for structural, all available P goes to structural */
    if (f->p_surf_struct_litter > psurf)
      f->p_surf_struct_litter = psurf;
    if (f->p_soil_struct_litter > psoil)
      f->p_soil_struct_litter = psoil;
  }

  /* remaining P goes to metabolic pools */
  f->p_surf_metab_litter = psurf - f->p_surf_struct_litter;
  f->p_soil_metab_litter = psoil - f->p_soil_struct_litter;

  return;
}

void pfluxes_from_structural_pools(fluxes *f, params *p, state *s) {
  /* from structural pool */
  double sigwt;
  double structout_surf = s->structsurfp * p->decayrate[0];
  double structout_soil = s->structsoilp * p->decayrate[2];

  sigwt = structout_surf / (p->ligshoot * 0.7 + (1.0 - p->ligshoot) * 0.55);

  /* P flux from surface structural pool -> slow pool */
  f->p_surf_struct_to_slow = sigwt * p->ligshoot * 0.7;

  /* P flux surface structural pool -> active pool */
  f->p_surf_struct_to_active = sigwt * (1.0 - p->ligshoot) * 0.55;

  sigwt = structout_soil / (p->ligroot * 0.7 + (1. - p->ligroot) * 0.45);


  /* P flux from soil structural pool -> slow pool */
  f->p_soil_struct_to_slow = sigwt * p->ligroot * 0.7;

  /* N flux from soil structural pool -> active pool */
  f->p_soil_struct_to_active = sigwt * (1.0 - p->ligroot) * 0.45;

  return;
}

void pfluxes_from_metabolic_pool(fluxes *f, params *p, state *s) {

  /* P flux surface metabolic pool -> active pool */
  f->p_surf_metab_to_active = s->metabsurfp * p->decayrate[1];

  /* P flux soil metabolic pool  -> active pool */
  f->p_soil_metab_to_active = s->metabsoilp * p->decayrate[3];

  return;
}


void pfluxes_from_active_pool(fluxes *f, params *p, state *s,
                              double frac_microb_resp) {

  double activeout, sigwt;
  /* P fluxes from active pool */
  activeout = s->activesoilp * p->decayrate[4];
  sigwt = activeout / (1.0 - frac_microb_resp);

  /* P flux active pool -> slow pool */
  f->p_active_to_slow = sigwt * (1.0 - frac_microb_resp - 0.004);

  /* P flux active pool -> passive pool */
  f->p_active_to_passive = sigwt * 0.004;

  return;
}

void pfluxes_from_slow_pool(fluxes *f, params *p, state *s) {
  /* P fluxes from slow pools */

  double slowout = s->slowsoilp * p->decayrate[5];
  double sigwt = slowout / 0.45;

  //fprintf(stderr, "slowsoilp %f\n", s->slowsoilp);

  /* P flux slow pool -> active pool */
  f->p_slow_to_active = sigwt * 0.42;

  /* slow pool -> passive pool */
  f->p_slow_to_passive = sigwt * 0.03;

  return;
}

void pfluxes_from_passive_pool(fluxes *f, params *p, state *s) {
  /* P fluxes from passive pool */

  /* P flux passive pool -> active pool */
  f->p_passive_to_active = s->passivesoilp * p->decayrate[6];

  return;
}

void calculate_p_parent_fluxes(fluxes *f, params *p, state *s) {
  /* Calculate weathering of parent P materials, i.e.
     the fluxes enterring into mineral P pool;

     Fluxes in = out so that parent P pool is a constant pool;

     Parameters:

     Return:
   p_atm_dep: atmospheric P deposition
   p_par_to_min: parent to mineral P pool
  */

  /* atmospheric P deposition rate */
  f->p_atm_dep = p->p_atm_deposition;

  //fprintf(stderr, "p_atm_deposition %f\n", p->p_atm_deposition);

  /* parent material weathering */
  f->p_par_to_min = p->p_rate_par_weather * s->inorgparp;

  //fprintf(stderr, "p_par_to_min %f\n", f->p_par_to_min);

  return;

}

void calculate_p_mineralisation(fluxes *f) {
  /* P gross mineralisation rate is given by the excess of P outflows
  over inflows. P mineralisation is the process by which organic P is
  converted to plant available inorganic P, i.e. the microbial conversion
  of organic P to H2PO4- or HPO42- forms of plant available P known as
  orthophosphates.

  Returns:
  --------
  value : float
  Gross P mineralisation
  Unit: t/ha/d
  */
  f->pgross =  (f->p_surf_struct_to_slow + f->p_surf_struct_to_active +
                f->p_soil_struct_to_slow + f->p_soil_struct_to_active +
                f->p_surf_metab_to_active + f->p_soil_metab_to_active +
                f->p_active_to_slow + f->p_active_to_passive +
                f->p_slow_to_active + f->p_slow_to_passive +
                f->p_passive_to_active);

  return;
}

void calculate_p_immobilisation(fluxes *f, params *p, state *s, double *pimmob,
                                double *active_pc_slope, double *slow_pc_slope,
                                double *passive_pc_slope) {
  /* P immobilised in new soil organic matter, the reverse of
  mineralisation. Micro-organisms in the soil compete with plants for P.
  Immobilisation occurs when plant available P forms are consumed by microbes, turning
  the P into organic P forms that are not available to plants.

  When C:P ratio is high the microorganisms need more P from
  the soil to decompose the carbon in organic materials. This P will be
  immobilised until these microorganisms die and the P is
  released.

  General equation for new soil P:C ratio vs Pmin, expressed as linear
  equation passing through point Pmin0, actpcmin (etc). Values can be
  Pmin0=0, Actpc0=Actpcmin

  if Pmin < Pmincrit:
  New soil P:C = soil P:C (when Pmin=0) + slope * Pmin

  if Pmin > Pmincrit
  New soil P:C = max soil P:C

  NB P:C ratio of new passive SOM can change even if assume Passiveconst

  Returns:
  --------
  pimob : float
  P immobilsed
  */
  double pmin, arg1, arg2, arg3, numer1, numer2, denom;

  /* P:C new SOM - active, slow and passive */
  *active_pc_slope = calculate_pc_slope(p, p->actpcmax, p->actpcmin);
  *slow_pc_slope = calculate_pc_slope(p, p->slowpcmax, p->slowpcmin);
  *passive_pc_slope = calculate_pc_slope(p, p->passpcmax, p->passpcmin);

  /* convert units */
  pmin = p->pmin0 / M2_AS_HA * G_AS_TONNES;

  arg1 = (p->passpcmin - *passive_pc_slope * pmin) * f->c_into_passive;
  arg2 = (p->slowpcmin - *slow_pc_slope * pmin) * f->c_into_slow;
  arg3 = f->c_into_active * (p->actpcmin - *active_pc_slope * pmin);
  numer1 = arg1 + arg2 + arg3;

  arg1 = f->c_into_passive * p->passpcmax;
  arg2 = f->c_into_slow * p->slowpcmax;
  arg3 = f->c_into_active * p->actpcmax;
  numer2 = arg1 + arg2 + arg3;

  arg1 = f->c_into_passive * *passive_pc_slope;
  arg2 = f->c_into_slow * *slow_pc_slope;
  arg3 = f->c_into_active * *active_pc_slope;
  denom = arg1 + arg2 + arg3;

  /* evaluate P immobilisation in new SOM */
  *pimmob = numer1 + denom * s->inorglabp;
  if (*pimmob > numer2)
    *pimmob = numer2;

  return;
}


void calc_p_net_mineralisation(fluxes *f) {
  /* P Net mineralisation from microbial activity,
  excluding the (- f->p_sorb_to_ssorb + f->p_ssorb_to_sorb activity) */
  f->pmineralisation = f->pgross - f->pimmob + f->plittrelease;

  //fprintf(stderr, "plitrelease %f\n", f->plittrelease);
  return;
}

double calculate_pc_slope(params *p, double pcmax, double pcmin) {
  /* Returns P:C ratio of the mineral pool slope
  Need to check back for good relationships; Currently using olde NC relationship
  based on Parton et al., 1993

  Parameters
  ----------
  pcmax : float
  SOM pools maximum P:C
  pcmin: float
  SOM pools minimum P:C

  Returns:
  --------
  value : float
  SOM pool P:C ratio
  */
  double arg1, arg2, conv;

  arg1 = pcmax - pcmin;
  arg2 = p->pmincrit - p->pmin0;
  conv = M2_AS_HA / G_AS_TONNES;

  return (arg1 / arg2 * conv);
}


void calculate_p_biochemical_mineralisation(fluxes *f, params *p,
                                            state *s) {
  /*
  Calculate biochemical P mineralisation rate for slow SOM pool only;

  Ref: Wang et al., 2007, GB1018.

  Parameters
  ----------
  n_cost_of_p: float
  N cost of plant root P uptake [g N (g P)-1]

  crit_n_cost_of_p: float
  The critical value of N cost of root P uptake above which
  phosphatase production starts [g N (g P)-1]

  max_p_biochemical: float
  Maximum rate of biochemical P mineralisation [g P m-2 y-1]

  biochemical_p_constant: float
  Michaelis-Menten constant for biochemical P mineralisation [g N (g P)-1]

  c_gain_of_p: float
  Carbon per unit of P uptake (NPP/Puptake) [g C (g P)-1]

  c_gain_of_n: float
  Carbon per unit of N uptake (NPP/Nuptake) [g C (g N)-1]

  Return
  ----------
  p_slow_biochemical: float
  biochemical P mineralisation rate in slow SOM pool
  Unit [t/ha];

  */

  double n_cost_of_p, c_gain_of_p, c_gain_of_n;

  /* Calculate C gain per unit P */
  if (f->puptake > 0.0) {
    c_gain_of_p = f->npp/f->puptake;
  } else {
    c_gain_of_p = 0.0;
  }

  /* Calculate C gain per unit N */
  if (f->nuptake > 0.0) {
    c_gain_of_n = f->npp/f->nuptake;
  } else {
    c_gain_of_n = 0.0;
  }

  /* Calculate n cost of p */
  if(c_gain_of_n > 0.0) {
    n_cost_of_p = c_gain_of_p / c_gain_of_n;
  } else {
    n_cost_of_p = 0.0;
  }

  /* Phosphatase production start when N cost of P is greater than critical
  N cost of P value, and when carbon uptake is more strongly limited by the uptake
  of P than by the uptake of N */
  if(c_gain_of_p > c_gain_of_n && n_cost_of_p > p->crit_n_cost_of_p) {
    f->p_slow_biochemical = p->max_p_biochemical *
                            ((n_cost_of_p - p->crit_n_cost_of_p) /
                            (n_cost_of_p - p->crit_n_cost_of_p + p->biochemical_p_constant));
  } else {
    f->p_slow_biochemical = 0.0;
  }

  return;
}

void calculate_p_min_fluxes(fluxes *f, params *p, state *s) {
  /* Calculate the mineral P fluxes (in and out)

  Parameters:
  --------

  Returns:
  --------
  value : float
  p_lab_in: P fluxes added into lab P;
  p_lab_out: P fluxes out lab P;
  p_sorb_in
  p_sorb_out

  */
  double numer, denom1, denom2;
  double min_frac_p_available_to_plant = 0.4;
  double max_frac_p_available_to_plant = 0.8;
  double mineral_n_with_max_p = 0.02;              /* Unit [t N ha-1] */
  double tot_in, tot_out;

  //fprintf(stderr, "flag 5 p_min_fluxes \n");

  /* Note: pmineralisation can be negative, and therefore, during spinup, when
  sorbp stock is low, the current approach resulted in negative sorbp in some cases,
  but it does not affect the final equilibrated state, so leave as is until better method
  is found */
  tot_in = f->p_par_to_min + f->pmineralisation +
           f->purine + f->p_slow_biochemical +
           f->p_ssorb_to_min;

  if(s->inorglabp > 0) {
    tot_out = f->puptake + f->ploss + f->p_min_to_ssorb;
  } else {
    f->puptake = 0.0;
    f->ploss = 0.0;
    f->p_min_to_ssorb = 0.0;
    tot_out = f->puptake + f->ploss + f->p_min_to_ssorb;
  }

  /* Use soil order to obtain smax and ks values */
  soil_soprtion_parameters(p->soil_order, p);

  /* Calculate lab P dynamics */
  numer = p->smax * p->ks;
  denom1 = (s->inorglabp + p->ks) * (s->inorglabp + p->ks);
  f->p_lab_in = tot_in / (1.0 + numer/denom1);
  f->p_lab_out = tot_out / (1.0 + numer/denom1);

  /* calculate sorb P dynamics */
  denom2 = (s->inorglabp + p->ks) * (s->inorglabp + p->ks) + numer;
  f->p_sorb_in = tot_in * (numer / denom2);
  f->p_sorb_out = tot_out * (numer / denom2);

  //fprintf(stderr, "mass balance p_min 1 %f\n",
  //        ((f->p_lab_in + f->p_sorb_in) - tot_in)*100000000000000);
  //fprintf(stderr, "mass balance p_min 2 %f\n",
  //        ((f->p_lab_out + f->p_sorb_out) - tot_out)*100000000000000);

  /* calculating fraction of labile P available for plant uptake */
  p->p_lab_avail = MAX(min_frac_p_available_to_plant,
                       MIN(min_frac_p_available_to_plant + s->inorgn *
                       (max_frac_p_available_to_plant - min_frac_p_available_to_plant) /
                       mineral_n_with_max_p, max_frac_p_available_to_plant));

  return;

}

void calculate_p_ssorb_to_sorb(state *s, fluxes *f, params *p, control *c) {
  /*calculate P transfer from strongly sorbed P pool to
   sorbed P pool;

   Parameters
   ----------
   phtextint: float
   intercept for the texture equation of strongly sorbed P depends upon
   pH input;

   phtextslope: float
   slope value used in determining effect of sand content on ssorb P flow
   to mineral P;

   psecmn: float
   controls the flow from secondary to mineral P;

   Returns:
   ----------
   p_ssorb_to_min: float
   flux rate of p strongly sorbed pool to p sorbed pool;

   */

  double phtextint;
  double dely, delx, xslope, yint;
  int cntrl_text_p = c->text_effect_p;

  if (cntrl_text_p == 1) {

    dely = p->phtextmax - p->phtextmin;
    delx = p->phmax - p->phmin;

    xslope = dely/delx;
    yint = p->phtextmin - xslope * p->phmin;

    if (p->soilph < p->phmin) {
      phtextint = p->phtextmin;
    } else if (p->soilph > p->phmax) {
      phtextint = p->phtextmax;
    } else {
      phtextint = xslope * p->soilph + yint;
    }

    f->p_ssorb_to_min = MAX(0.0, (phtextint + p->phtextslope * (1.0 - p->finesoil)) *
                         s->inorgssorbp);

  } else {
    if (s->inorgssorbp > 0.0) {
      f->p_ssorb_to_min = p->psecmnp * s->inorgssorbp;
    } else {
      f->p_ssorb_to_min = 0.0;
    }

  }

  return;
}

void calculate_p_sorb_to_ssorb(state *s, fluxes *f, params *p) {

  /* P flux from sorbed pool to strongly sorbed P pool */
  if (s->inorgsorbp > 0.0) {
  f->p_min_to_ssorb = p->rate_sorb_ssorb * s->inorgsorbp;
  } else {
    f->p_min_to_ssorb = 0.0;
  }

  return;
}

void calculate_p_ssorb_to_occ(state *s, fluxes *f, params *p) {

  /* P flux from strongly sorbed pool to occluded P pool */
  if (s->inorgssorbp > 0.0) {
  f->p_ssorb_to_occ = p->rate_ssorb_occ * s->inorgssorbp;
  } else {
    f->p_ssorb_to_occ = 0.0;
  }

  return;
}


void calculate_ppools(control *c, fluxes *f, params *p, state *s,
                      double active_pc_slope, double slow_pc_slope,
                      double passive_pc_slope) {
  /*
  Update P pools in the soil

  Parameters
  ----------
  active_pc_slope : float
  active PC slope
  slow_pc_slope: float
  slow PC slope
  passive_pc_slope : float
  passive PC slope

  */

  double p_into_active, p_out_of_active, p_into_slow, p_out_of_slow,
  p_into_passive, p_out_of_passive, arg, active_pc, fixp, slow_pc,
  pass_pc;

  double net_parent;

  /*
  net P release implied by separation of litter into structural
  & metabolic. The following pools only fix or release P at their
  limiting p:c values.
  */

  /* P released or fixed from the P inorganic labile pool is incremented with
  each call to pc_limit and stored in f->plittrelease */
  f->plittrelease = 0.0;

  s->structsurfp += (f->p_surf_struct_litter -
    (f->p_surf_struct_to_slow +
    f->p_surf_struct_to_active));

  s->structsoilp += (f->p_soil_struct_litter -
    (f->p_soil_struct_to_slow + f->p_soil_struct_to_active));

  if (c->strpfloat == 0) {
    s->structsurfp += pc_limit(f, s->structsurf, s->structsurfp,
                               1.0/p->structcp, 1.0/p->structcp);
    s->structsoilp += pc_limit(f, s->structsoil, s->structsoilp,
                               1.0/p->structcp, 1.0/p->structcp);
  }

  s->metabsurfp += f->p_surf_metab_litter - f->p_surf_metab_to_active;
  s->metabsurfp += pc_limit(f, s->metabsurf, s->metabsurfp,
                            1.0/150.0, 1.0/80.0);                                  /* pcmin & pcmax from Parton 1989 fig 2 */

  s->metabsoilp += (f->p_soil_metab_litter - f->p_soil_metab_to_active);
  s->metabsoilp += pc_limit(f, s->metabsoil, s->metabsoilp,
                            1.0/150.0, 1.0/80.0);                                  /* pcmin & pcmax from Parton 1989 fig 2 */


  /* When nothing is being added to the metabolic pools, there is the
  potential scenario with the way the model works for tiny bits to be
  removed with each timestep. Effectively with time this value which is
  zero can end up becoming zero but to a silly decimal place */
  precision_control_soil_p(f, s);

  /* Update SOM pools */
  p_into_active = (f->p_surf_struct_to_active + f->p_soil_struct_to_active +
  f->p_surf_metab_to_active + f->p_soil_metab_to_active +
  f->p_slow_to_active + f->p_passive_to_active);

  p_out_of_active = f->p_active_to_slow + f->p_active_to_passive;

  p_into_slow = (f->p_surf_struct_to_slow + f->p_soil_struct_to_slow +
    f->p_active_to_slow);

  p_out_of_slow = f->p_slow_to_active + f->p_slow_to_passive + f->p_slow_biochemical;
  p_into_passive = f->p_active_to_passive + f->p_slow_to_passive;
  p_out_of_passive = f->p_passive_to_active;

  /* P:C of the SOM pools increases linearly btw prescribed min and max
  values as the Pconc of the soil increases. */
  arg = s->inorglabp - p->pmin0 / M2_AS_HA * G_AS_TONNES;

  /* active */
  active_pc = p->actpcmin + active_pc_slope * arg;
  if (active_pc > p->actpcmax)
    active_pc = p->actpcmax;

  /* release P to Inorganic labile pool or fix P from the Inorganic pool in order
  to normalise the P:C ratio of a net flux */
  fixp = pc_flux(f->c_into_active, p_into_active, active_pc);
  s->activesoilp += p_into_active + fixp - p_out_of_active;

  /* slow */
  slow_pc = p->slowpcmin + slow_pc_slope * arg;
  if (slow_pc > p->slowpcmax)
    slow_pc = p->slowpcmax;

  /* release P to Inorganic pool or fix P from the Inorganic pool in order
  to normalise the P:C ratio of a net flux */
  fixp = pc_flux(f->c_into_slow, p_into_slow, slow_pc);
  s->slowsoilp += p_into_slow + fixp - p_out_of_slow;

  /* passive, update passive pool only if passiveconst=0 */
  pass_pc = p->passpcmin + passive_pc_slope * arg;
  if (pass_pc > p->passpcmax)
    pass_pc = p->passpcmax;

  /* release P to Inorganic pool or fix P from the Inorganic pool in order
  to normalise the P:C ratio of a net flux */
  fixp = pc_flux(f->c_into_passive, p_into_passive, pass_pc);
  s->passivesoilp += p_into_passive + fixp - p_out_of_passive;

  //fprintf(stderr, "inorglabp 1 %f\n", s->inorglabp);

  /* Daily increment of soil inorganic labile and sorbed P pool */
  s->inorglabp += f->p_lab_in - f->p_lab_out;
  s->inorgsorbp += f->p_sorb_in - f->p_sorb_out;

  //fprintf(stderr, "psorb calc %f\n", (9 * s->inorglabp)/(0.0012+s->inorglabp));

  //fprintf(stderr, "inorglabp %f\n", s->inorglabp);
  //fprintf(stderr, "inorgsorbp %f\n", s->inorgsorbp);

  /* Daily increment of soil inorganic available P pool (lab + sorb) */
  s->inorgavlp = s->inorglabp + s->inorgsorbp;

  /* Daily increment of soil inorganic secondary P pool (strongly sorbed) */
  s->inorgssorbp += f->p_min_to_ssorb - f->p_ssorb_to_occ - f->p_ssorb_to_min;


  /* Daily increment of soil inorganic occluded P pool */
  s->inorgoccp += f->p_ssorb_to_occ;

  /* Daily increment of soil inorganic parent P pool */
  s->inorgparp += f->p_atm_dep - f->p_par_to_min;

  //fprintf(stderr, "flag 6 ppools \n");

  return;
}


double pc_limit(fluxes *f, double cpool, double ppool, double pcmin,
                double pcmax) {
  /* Release P to 'Inorglabp' pool or fix P from 'Inorglabp', in order to keep
  the  P:C ratio of a litter pool within the range 'pcmin' to 'pcmax'.

  Parameters:
  -----------
  cpool : float
  various C pool (state)
  ppool : float
  various P pool (state)
  pcmin : float
  min P:C ratio
  pcmax : float
  max P:C ratio

  Returns:
  --------
  fix/rel : float
  amount of P to be added/released from the inorganic pool

  */
  double rel, fix;
  double pmax = cpool * pcmax;
  double pmin = cpool * pcmin;

  //fprintf(stderr, "ppool %f\n", ppool*100000);
  //fprintf(stderr, "pmax %f\n", pmax*100000);
  //fprintf(stderr, "pmin %f\n", pmin*100000);

  if (ppool > pmax) {
    /* release */
    rel = ppool - pmax;
    f->plittrelease += rel;
    return (-rel);
  } else if (ppool < pmin) {
    /* fix */
    fix = pmin - ppool;
    f->plittrelease -= fix;
    return (fix);
  } else {
    return (0.0);
  }
}

double pc_flux(double cflux, double pflux, double pc_ratio) {
  /*
  Release P to Inorganic pool or fix P from the Inorganic pool in order
  to normalise the P:C ratio of a net flux

  Parameters:
  -----------
  cflux : float
  C flux into SOM pool
  pflux : float
  P flux into SOM pool
  pc_ratio : float
  preferred P:C ratio

  Returns:
  fix : float
  Returns the amount of P required to be fixed
  */

  return (cflux * pc_ratio) - pflux;
}


void precision_control_soil_p(fluxes *f, state *s) {
  /* Detect very low values in state variables and force to zero to
  avoid rounding and overflow errors */

  double tolerance = 1E-08, excess;

  if (s->metabsurfp < tolerance) {
    excess = s->metabsurfp;
    f->p_surf_metab_to_active = excess;
    s->metabsurfp = 0.0;
  }

  if (s->metabsoilp < tolerance) {
    excess = s->metabsoilp;
    f->p_soil_metab_to_active = excess;
    s->metabsoilp = 0.0;
  }

  return;
}
