/* ============================================================================
* Soil C and N flows into 4 litter pools (structural and metabolic, both
* above and belowground) and 3 SOM pools (Active, slow and passive). In
* essence the CENTURY model.
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
*   22.02.2015
*
* =========================================================================== */
#include "soils.h"

void calculate_csoil_flows(control *c, fluxes *f, params *p, state *s,
                           double tsoil) {
    double lnleaf, lnroot, nc_leaf_litter;
    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);

    soil_temp_factor(f, tsoil);

    /* calculate model decay rates */
    calculate_decay_rates(f, p, s);

    /*
       plant litter inputs to the metabolic and structural pools determined
       by ratio of lignin/N ratio
    */
    ligin_nratio(c, f, p, &lnleaf, &lnroot);
    p->fmleaf = metafract(lnleaf);
    p->fmroot = metafract(lnroot);
    nc_leaf_litter = ratio_of_litternc_to_live_leafnc(c, f, p);

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
    adfac = s->wtfac_topsoil * f->tfac_soil_decomp;

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

void soil_temp_factor(fluxes *f, double tsoil) {
    /*
    Soil-temperature activity factor (A9). Fit to Parton's fig 2a

    Parameters:
    -----------
    tsoil : double
        soil temperature (deg C)

    Returns:
    --------
    tfac : float
        soil temperature factor [degC]

    */

    if (tsoil > 0.0) {
        f->tfac_soil_decomp = (0.0326 + 0.00351 * pow(tsoil, 1.652) -
                                pow((tsoil / 41.748), 7.19));
        if (f->tfac_soil_decomp < 0.0)
        f->tfac_soil_decomp = 0.0;
    } else {
        /* negative number cannot be raised to a fractional power
           number would need to be complex */
        f->tfac_soil_decomp = 0.0;
    }

    return;
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


void ligin_nratio(control *c, fluxes *f, params *p, double *lnleaf,
                  double *lnroot) {
    /* Estimate Lignin/N ratio, as this dictates the how plant litter is
    seperated between metabolic and structural pools.

    Returns:
    --------
    lnleaf : float
        lignin:N ratio of leaf
    lnroot : float
        lignin:N ratio of fine root
    */
    double nc_leaf_litter, nc_root_litter;

    nc_leaf_litter = ratio_of_litternc_to_live_leafnc(c, f, p);
    nc_root_litter = ratio_of_litternc_to_live_rootnc(c, f, p);

    if (float_eq(nc_leaf_litter, 0.0))
        /* catch divide by zero if we have no leaves */
        *lnleaf = 0.0;
    else
        *lnleaf = p->ligshoot / p->cfracts / nc_leaf_litter;

    if (float_eq(nc_root_litter, 0.0))
        /* catch divide by zero if we have no roots */
        *lnroot = 0.0;
    else
        *lnroot = p->ligroot / p->cfracts / nc_root_litter;

    return;

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

    /*
    ** Surface (leaves, branches, stem) Litter
    */

    /* ...to the structural pool*/
    f->surf_struct_litter = (f->deadleaves * (1.0 - p->fmleaf) +
                             f->deadbranch + f->deadstems + f->faecesc *
                             (1.0 - p->fmfaeces));

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

    /* C fluxes from structural pools */

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
    /* C fluxes from metabolic pools */

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
    /* C fluxes from active pools */

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
    /* C fluxes from slow pools */

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

    /* store the C SOM fluxes for Nitrogen calculations */
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
                           double ndep) {

    /* Fraction of C lost due to microbial respiration */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double nsurf, nsoil, active_nc_slope, slow_nc_slope, passive_nc_slope;

    grazer_inputs(c, f, p);
    inputs_from_plant_litter(f, p, &nsurf, &nsoil);
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

    /* Update model soil N pools */
    calculate_npools(c, f, p, s, active_nc_slope, slow_nc_slope,
                     passive_nc_slope, ndep);

    /* calculate N net mineralisation */
    calc_net_mineralisation(f);

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

void inputs_from_plant_litter(fluxes *f, params *p, double *nsurf,
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
    if (c->strfloat) {

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


void calc_net_mineralisation(fluxes *f) {
    /* N Net mineralisation from microbial activity
       i.e. excess of N outflows over inflows */
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
                      double passive_nc_slope, double ndep) {
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

    if (c->strfloat == FALSE) {
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
    s->inorgn += ((f->ngross + ndep + f->nurine - f->nimmob -
                   f->nloss - f->nuptake) + f->nlittrelease);

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
