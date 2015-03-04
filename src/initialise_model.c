#include "initialise_model.h"


void initialise_control(control *c) {
    /*
    *** Default values for control structure.
    */

    /* Values set via param file */
    strcpy(c->git_hash, "Err");

    c->ifp = NULL;
    c->ofp = NULL;
    c->ofp_hdr = NULL;
    strcpy(c->cfg_fname, "*NOT SET*");
    strcpy(c->met_fname, "*NOT SET*");
    strcpy(c->out_fname, "*NOT SET*");
    strcpy(c->out_fname_hdr, "*NOT SET*");
    strcpy(c->out_param_fname, "*NOT SET*");

    c->alloc_model = ALLOMETRIC;    /* C allocation scheme: FIXED, GRASSES, ALLOMETRIC */
    c->assim_model = MATE;          /* Photosynthesis model: BEWDY (not coded :p) or MATE */
    c->calc_sw_params = FALSE;      /* false=user supplies field capacity and wilting point, true=calculate them based on cosby et al. */
    c->deciduous_model = FALSE;     /* evergreen_model=False, deciduous_model=True */
    c->fixed_stem_nc = TRUE;        /* False=vary stem N:C with foliage, True=fixed stem N:C */
    c->fixleafnc = FALSE;           /* fixed leaf N C ? */
    c->grazing = FALSE;             /* Is foliage grazed? 0=No, 1=daily, 2=annual and then set disturbance_doy=doy */
    c->gs_model = MEDLYN;           /* Stomatal conductance model, currently only this one is implemented */
    c->model_optroot = FALSE;       /* Ross's optimal root model...not sure if this works yet...0=off, 1=on */
    c->modeljm = 2;                 /* modeljm=0, Jmax and Vcmax parameters are read in, modeljm=1, parameters are calculated from leaf N content, modeljm=2, Vcmax is calculated from leaf N content but Jmax is related to Vcmax */
    c->ncycle = TRUE;               /* Nitrogen cycle on or off? */
    c->nuptake_model = 2;           /* 0=constant uptake, 1=func of N inorgn, 2=depends on rate of soil N availability */
    c->output_ascii = TRUE;         /* If this is false you get a binary file as an output. */
    c->passiveconst = FALSE;        /* hold passive pool at passivesoil */
    c->print_options = DAILY;       /* DAILY=every timestep, END=end of run */
    c->ps_pathway = C3;             /* Photosynthetic pathway, c3/c4 */
    c->respiration_model = FIXED;   /* Plant respiration ... Fixed, TEMPERATURE or BIOMASS */
    c->strfloat = 0;                /* Structural pool input N:C varies=1, fixed=0 */
    c->sw_stress_model = 1;         /* JULES type linear stress func, or Landsberg and Waring non-linear func */
    c->trans_model = 0;             /* 0=trans from WUE, 1=Penman-Monteith, 2=Priestley-Taylor */
    c->use_eff_nc = 0;              /* use constant leaf n:c for  metfrac s */
    c->water_stress = TRUE;         /* water stress modifier turned on=TRUE (default)...ability to turn off to test things without drought stress = FALSE */
    c->spin_up = FALSE;             /* Spin up to a steady state? If False it just runs the model */

    /* Internal calculated */
    c->num_years = 0;               /* Total number of years simulated */
    c->num_days = 0;                /* Number of days in a year: 365/366 */
    c->PRINT_GIT = FALSE;           /* print the git hash to the cmd line and exit? Called from cmd line parsar */

    return;
}

void initialise_params(params *p) {
    /*
    *** Default values for params structure.
    */
    int i;
    p->actncmax = 0.333333;
    p->actncmin = 0.066667;
    p->ageold = 10000.0;
    p->ageyoung = 0.0;
    p->albedo = 0.123;
    p->alpha_c4 = 0.06;
    p->alpha_j = 0.26;
    p->b_root = -999.9;
    p->b_topsoil =-999.9;
    p->bdecay = 0.02;
    p->branch0 = 5.61;
    p->branch1 = 0.346;
    p->bretrans = 0.0;
    p->c_alloc_bmax = 0.1;
    p->c_alloc_bmin = 0.1;
    p->c_alloc_cmax = 0.0;
    p->c_alloc_fmax = 0.35;
    p->c_alloc_fmin = 0.15;
    p->c_alloc_rmax = 0.35;
    p->c_alloc_rmin = 0.05;
    p->cfracts = 0.5;
    p->crdecay = 0.0;
    p->cretrans = 0.0;
    p->croot0 = 0.34;
    p->croot1 = 0.84;
    p->ctheta_root = 0.4;
    p->ctheta_topsoil = 0.5;
    p->cue = 0.5;
    p->d0 = 0.0;
    p->d0x = 0.35;
    p->d1 = 0.0;
    p->delsj = 644.4338;
    p->density = 420.0;
    p->direct_frac = 0.5;
    p->displace_ratio = 0.78;
    p->disturbance_doy = 1.0;
    p->dz0v_dh = 0.075;
    p->eac = 79430.0;
    p->eag = 37830.0;
    p->eaj = 43790.0;
    p->eao = 36380.0;
    p->eav = 51560.0;
    p->edj = 200000.0;
    p->faecescn = 25.0;
    p->faecesn = 0.0;
    p->fdecay = 0.59988;
    p->fdecaydry = 0.59988;
    p->fhw = 0.8;
    p->finesoil = 0.51;
    p->fracfaeces = 0.3;
    p->fracteaten = 0.5;
    p->fractosoil = 0.85;
    p->fractup_soil = 0.5;
    p->fretrans = 0.5;
    p->g1 = 2.74;
    p->gamstar25 = 42.75;
    p->growth_efficiency = 0.7;
    p->height0 = 5.0;
    p->height1 = 30.0;
    p->heighto = 4.826;
    p->htpower = 0.35;
    p->intercep_frac = 0.15;
    p->jmax = -999.9;
    p->jmaxna = 62.0;
    p->jmaxnb = 0.0;
    p->jv_intercept = 0.0;
    p->jv_slope = 1.86;
    p->kc25 = 404.9;
    p->kdec1 = 3.965571;
    p->kdec2 = 14.61;
    p->kdec3 = 4.904786;
    p->kdec4 = 18.262499;
    p->kdec5 = 7.305;
    p->kdec6 = 0.198279;
    p->kdec7 = 0.006783;
    p->kext = 0.5;
    p->knl = 0.01;
    p->ko25 = 278400.0;
    p->kq10 = 0.08;
    p->kr = 0.5;
    p->lai_closed = 0.5;
    p->latitude = 35.9;
    p->leafsap0 = 8000.0;
    p->leafsap1 = 3060.0;
    p->ligfaeces = 0.25;
    p->ligroot = 0.22;
    p->ligshoot = 0.24;
    p->liteffnc = 0.0;
    p->max_intercep_lai = 3.0;
    p->measurement_temp = 25.0;
    p->ncbnew = 0.003;
    p->ncbnewz = 0.003;
    p->nccnew = 0.003;
    p->nccnewz = 0.003;
    p->ncmaxfold = 0.04;
    p->ncmaxfyoung = 0.04;
    p->ncmaxr = 0.03;
    p->ncrfac = 0.8;
    p->ncwimm = 0.003;
    p->ncwimmz = 0.003;
    p->ncwnew = 0.003;
    p->ncwnewz = 0.003;
    p->nf_crit = 0.015;
    p->nf_min = 0.005;
    p->nmax = 0.24;
    p->nmin = 0.95;
    p->nmin0 = 0.0;
    p->nmincrit = 2.0;
    p->ntheta_root = 3.0;
    p->ntheta_topsoil = 5.0;
    p->nuptakez = 0.0;
    p->oi = 205000.0;
    p->passivesoilnz = 1.0;
    p->passivesoilz = 1.0;
    p->passncmax = 0.142857;
    p->passncmin = 0.1;
    p->prescribed_leaf_NC = 0.03;
    p->previous_ncd = 35.0;
    p->psi_sat_root = -999.9;
    p->psi_sat_topsoil = -999.9;
    p->qs = 1.0; /* exponent in water stress modifier, =1.0 JULES type representation, the smaller the values the more curved the depletion. */
    p->r0 = 0.1325;
    p->rateloss = 0.5;
    p->rateuptake = 2.7;
    p->rdecay = 0.33333;
    p->rdecaydry = 0.33333;
    p->retransmob = 0.0;
    p->rfmult = 1.0;
    p->root_exu_CUE = -999.9;
    p->rooting_depth = 750.0;
    strcpy(p->rootsoil_type, "clay");
    p->rretrans = 0.0;
    p->sapturnover = 0.1;
    p->sla = 4.4;
    p->slamax = 4.4;
    p->slazero = 4.4;
    p->slowncmax = 0.066666;
    p->slowncmin = 0.025;
    p->store_transfer_len = -999.9;
    p->structcn = 150.0;
    p->structrat = 0.0;
    p->targ_sens = 0.5;
    p->theta = 0.7;
    p->theta_sat_root = -999.9;
    p->theta_sat_topsoil = -999.9;
    p->topsoil_depth = 350.0;
    strcpy(p->topsoil_type, "clay_loam");
    p->vcmax = -999.9;
    p->vcmaxna = 22.29;
    p->vcmaxnb = 8.45;
    p->watdecaydry = 0.0;
    p->watdecaywet = 0.1;
    p->wcapac_root = 96.75;
    p->wcapac_topsoil = 25.8;
    p->wdecay = 0.02;
    p->wetloss = 0.5;
    p->wretrans = 0.0;
    p->z0h_z0m = 1.0;
    p->fmleaf = 0.0;
    p->fmroot = 0.0;
    p->fmfaeces = 0.0;
    p->growing_seas_len = 0;

    for (i = 0; i < 7; i++) {
        p->decayrate[i] = 0.0;
    }
}


void initialise_fluxes(fluxes *f) {
    /*
    ** Default values for fluxes structure.
    */
    int i = 0;

    /* C fluxes */
    f->gpp_gCm2 = 0.0;
    f->npp_gCm2 = 0.0;
    f->gpp = 0.0;
    f->npp = 0.0;
    f->nep = 0.0;
    f->auto_resp = 0.0;
    f->hetero_resp = 0.0;
    f->retrans = 0.0;
    f->apar = 0.0;

    /* N fluxes */
    f->nuptake = 0.0;
    f->nloss = 0.0;
    f->npassive = 0.0;              /* n passive -> active */
    f->ngross = 0.0;                /* N gross mineralisation */
    f->nimmob = 0.0;                /* N immobilisation in SOM */
    f->nlittrelease = 0.0;          /* N rel litter = struct + metab */
    f->activelossf = 0.0;           /* frac of active C -> CO2 */
    f->nmineralisation = 0.0;

    /* water fluxes */
    f->wue = 0.0;
    f->et = 0.0;
    f->soil_evap = 0.0;
    f->transpiration = 0.0;
    f->erain = 0.0;
    f->interception = 0.0;
    f->runoff = 0.0;
    f->gs_mol_m2_sec = 0.0;
    f->ga_mol_m2_sec = 0.0;
    f->omega = 0.0;

    /* daily C production */
    f->cpleaf = 0.0;
    f->cproot = 0.0;
    f->cpcroot = 0.0;
    f->cpbranch = 0.0;
    f->cpstem = 0.0;

    /* daily N production */
    f->npleaf = 0.0;
    f->nproot = 0.0;
    f->npcroot = 0.0;
    f->npbranch = 0.0;
    f->npstemimm = 0.0;
    f->npstemmob = 0.0;

    /* dying stuff */
    f->deadleaves = 0.0;   /* Leaf litter C production (t/ha/yr) */
    f->deadroots = 0.0;    /* Root litter C production (t/ha/yr) */
    f->deadcroots = 0.0;   /* Coarse root litter C production (t/ha/yr) */
    f->deadbranch = 0.0;   /* Branch litter C production (t/ha/yr) */
    f->deadstems = 0.0;    /* Stem litter C production (t/ha/yr) */
    f->deadleafn = 0.0;    /* Leaf litter N production (t/ha/yr) */
    f->deadrootn = 0.0;    /* Root litter N production (t/ha/yr) */
    f->deadcrootn = 0.0;   /* Coarse root litter N production (t/ha/yr) */
    f->deadbranchn = 0.0;  /* Branch litter N production (t/ha/yr) */
    f->deadstemn = 0.0;    /* Stem litter N production (t/ha/yr) */
    f->deadsapwood = 0.0;

    /* grazing stuff */
    f->ceaten = 0.0;       /* C consumed by grazers (t C/ha/y) */
    f->neaten = 0.0;       /* N consumed by grazers (t C/ha/y) */
    f->faecesc = 0.0;      /* Flux determined by faeces C:N */
    f->nurine = 0.0;       /* Rate of N input to soil in urine (t/ha/y) */

    f->leafretransn = 0.0;

    /* C&N Surface litter */
    f->surf_struct_litter = 0.0;
    f->surf_metab_litter = 0.0;
    f->n_surf_struct_litter = 0.0;
    f->n_surf_metab_litter = 0.0;

    /* C&N Root Litter */
    f->soil_struct_litter = 0.0;
    f->soil_metab_litter = 0.0;
    f->n_soil_struct_litter = 0.0;
    f->n_soil_metab_litter = 0.0;

    /* C&N litter fluxes to slow pool */
    f->surf_struct_to_slow = 0.0;
    f->soil_struct_to_slow = 0.0;
    f->n_surf_struct_to_slow = 0.0;
    f->n_soil_struct_to_slow = 0.0;

    /* C&N litter fluxes to active pool */
    f->surf_struct_to_active = 0.0;
    f->soil_struct_to_active = 0.0;
    f->n_surf_struct_to_active = 0.0;
    f->n_soil_struct_to_active = 0.0;

    /* Metabolic fluxes to active pool */
    f->surf_metab_to_active = 0.0;
    f->soil_metab_to_active = 0.0;
    f->n_surf_metab_to_active = 0.0;
    f->n_soil_metab_to_active = 0.0;

    /* C fluxes out of active pool */
    f->active_to_slow = 0.0;
    f->active_to_passive = 0.0;
    f->n_active_to_slow = 0.0;
    f->n_active_to_passive = 0.0;

    /* C&N fluxes from slow to active pool */
    f->slow_to_active = 0.0;
    f->slow_to_passive = 0.0;
    f->n_slow_to_active = 0.0;
    f->n_slow_to_passive = 0.0;

    /* C&N fluxes from passive to active pool */
    f->passive_to_active = 0.0;
    f->n_passive_to_active = 0.0;

    /* C & N source fluxes from the active, slow and passive pools */
    f->c_into_active = 0.0;
    f->c_into_slow = 0.0;
    f->c_into_passive = 0.0;

    /* CO2 flows to the air */
    /* C flows to the air */
    for (i = 0; i < 7; i++) {
        f->co2_to_air[i] = 0.0;
    }

    /* C allocated fracs  */
    f->alleaf = 0.0;
    f->alroot = 0.0;
    f->alcroot = 0.0;
    f->albranch = 0.0;
    f->alstem = 0.0;

    /* Misc stuff */
    f-> cica_avg = 0.0; /* used in water balance, only when running mate model */

    f->rabove = 0.0;
    f->tfac_soil_decomp = 0.0;
    f->co2_rel_from_surf_struct_litter = 0.0;
    f->co2_rel_from_soil_struct_litter = 0.0;
    f->co2_rel_from_surf_metab_litter = 0.0;
    f->co2_rel_from_soil_metab_litter = 0.0;
    f->co2_rel_from_active_pool = 0.0;
    f->co2_rel_from_slow_pool = 0.0;
    f->co2_rel_from_passive_pool = 0.0;

    return;
}



void initialise_state(state *s) {

    /*
    *** Default values for state structure.
    */

    s->activesoil = 2.53010543182;
    s->activesoiln = 0.833516379296;
    s->age = 12.0;
    s->avg_albranch = 0.0;
    s->avg_alcroot = 0.0;
    s->avg_alleaf = 0.0;
    s->avg_alroot = 0.0;
    s->avg_alstem = 0.0;
    s->branch = 14.5137000708;
    s->branchn = 0.0442890661217;
    s->canht = 23.0964973582;
    s->croot = 0.0;
    s->crootn = 0.0;
    s->cstore = 0.01;
    s->inorgn = 0.0274523714275;
    s->max_lai = -999.9;
    s->max_shoot = -999.9;
    s->metabsoil = 0.135656771805;
    s->metabsoiln = 0.00542627087221;
    s->metabsurf = 0.0336324759951;
    s->metabsurfn = 0.0013452990398;
    s->nstore = 0.01;
    s->passivesoil = 59.5304597863;
    s->passivesoiln = 8.0134056319;
    s->pawater_root = 94.0415606424;
    s->pawater_topsoil = 24.7780118747;
    s->prev_sma = -999.9;
    s->root = 3.92887790342;
    s->root_depth = -9999.9;
    s->rootn = 0.076296932914;
    s->sapwood = 51.2600270003;
    s->shoot = 4.37991243755;
    s->shootn = 0.0978837857406;
    s->sla = 4.4;
    s->slowsoil = 46.8769593608;
    s->slowsoiln = 2.90664959452;
    s->stem = 87.6580936643;
    s->stemn = 0.263722246902;
    s->stemnimm = 0.263336697464;
    s->stemnmob = 0.00038554943772;
    s->structsoil = 0.917128200367;
    s->structsoiln = 0.00611418800245;
    s->structsurf = 7.10566198821;
    s->structsurfn = 0.0473710799214;

    s->wtfac_root = 1.0;
    return;
}
