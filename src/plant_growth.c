/* ============================================================================
* Calls photosynthesis model, water balance and evolves aboveground plant
* C & Nstate. Pools recieve C through allocation of accumulated photosynthate
* and N from both soil uptake and retranslocation within the plant. Key feedback
* through soil N mineralisation and plant N uptake
*
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   17.02.2015
*
* =========================================================================== */
#include "plant_growth.h"
#include "water_balance.h"


void calc_day_growth(canopy_wk *cw, control *c, fluxes *f, met_arrays *ma,
                     met *m, params *p, state *s, double day_length, int doy,
                     double fdecay, double rdecay)
{
   double previous_topsoil_store, dummy=0.0, previous_rootzone_store;
   double nitfac, pitfac, npitfac;
   double ncbnew, nccnew, ncwimm, ncwnew;
   double pcbnew, pccnew, pcwimm, pcwnew;
   int    recalc_wb;

    /* Store the previous days soil water store */
    previous_topsoil_store = s->pawater_topsoil;
    previous_rootzone_store = s->pawater_root;

    if (c->sub_daily) {
        /* calculate 30 min two-leaf GPP/NPP, respiration and water fluxes */
        canopy(cw, c, f, ma, m, p, s);
    } else {
        /* calculate daily GPP/NPP, respiration and update water balance */
        carbon_daily_production(c, f, m, p, s, day_length);
        calculate_water_balance(c, f, m, p, s, day_length, dummy, dummy, dummy);
    }

    /* leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
       foliage in young stand, and leaf P:C as a fraction of Pcmaxyoung;
    */
    nitfac = MIN(1.0, s->shootnc / p->ncmaxfyoung);
    pitfac = MIN(1.0, s->shootpc / p->pcmaxfyoung);
   
    /* checking for pcycle control parameter */ 
    if(c->pcycle == TRUE) {
       npitfac = MIN(nitfac, pitfac);
    } else {
       npitfac = nitfac;
    }

    /* figure out the C allocation fractions */
    if (c->deciduous_model){
        /* Allocation is annually for deciduous "tree" model, but we need to
           keep a check on stresses during the growing season and the LAI
           figure out limitations during leaf growth period. This also
           applies for deciduous grasses, need to do the growth stress
           calc for grasses here too. */
        if (s->leaf_out_days[doy] > 0.0) {

            calc_carbon_allocation_fracs(c, f, p, s, npitfac);

            /* store the days allocation fraction, we average these at the
               end of the year (for the growing season) */
            s->avg_alleaf += f->alleaf;
            s->avg_albranch += f->albranch;
            s->avg_alstem += f->alstem;
            s->avg_alroot += f->alroot;
            s->avg_alcroot += f->alcroot;

        }
    } else {
        /* daily allocation...*/
        calc_carbon_allocation_fracs(c, f, p, s, npitfac);
    }

    /* Distribute new C, N and P through the system */
    carbon_allocation(c, f, p, s, npitfac, doy);

    calculate_cnp_wood_ratios(c, p, s, npitfac, nitfac, pitfac, 
                              &ncbnew, &nccnew, &ncwimm,
                              &ncwnew, &pcbnew, &pccnew, &pcwimm,
                              &pcwnew);

    recalc_wb = np_allocation(c, f, p, s, ncbnew, nccnew, ncwimm, ncwnew,
                              pcbnew, pccnew, pcwimm, pcwnew,
                              fdecay, rdecay, doy);

    if (c->exudation && c->alloc_model != GRASSES) {
        calc_root_exudation(c, f, p, s);
    }

    /* If we didn't have enough N available to satisfy wood demand, NPP
       is down-regulated and thus so is GPP. We also need to recalculate the
       water balance given the lower GPP. */
    if (recalc_wb) {
        s->pawater_topsoil = previous_topsoil_store;
        s->pawater_root = previous_rootzone_store;

        if (c->sub_daily) {
            /* reduce transpiration to match cut back GPP
                -there isn't an obvious way to make this work at the 30 min
                 timestep, so invert T from WUE assumption and use that
                 to recalculate the end day water balance
            */
            f->transpiration = f->gpp_gCm2 / f->wue;
            update_water_storage_recalwb(c, f, p, s, m);

        } else {
            calculate_water_balance(c, f, m, p, s, day_length, dummy, dummy,
                                    dummy);
        }

    }
    update_plant_state(c, f, p, s, fdecay, rdecay, doy);
    precision_control(f, s);

    return;
}

void calc_root_exudation(control *c, fluxes *f, params *p, state *s) {
    /*
        Rhizodeposition (f->root_exc) is assumed to be a fraction of the
        current root growth rate (f->cproot), which increases with increasing
        N and P stress of the plant.
    */
    double CN_leaf, CN_ref, arg1;
    double CP_leaf, CP_ref, arg2;
    double frac_to_rexc;

    if (float_eq(s->shoot, 0.0) || float_eq(s->shootn, 0.0 || float_eq(s->shootp, 0.0))) {
        /* nothing happens during leaf off period */
        CN_leaf = 0.0;
        CP_leaf = 0.0;
        frac_to_rexc = 0.0;
    } else {

        if (c->deciduous_model) {
            /* broadleaf */
            CN_ref = 25.0;
            CP_ref = 600;   /* Needs to get empirical evidence for this number */
        } else {
            /* conifer */
            CN_ref = 42.0;
            CP_ref = 2500;  /* Needs to get empirical evidence for this number */
        }

        /*
        ** The fraction of growth allocated to rhizodeposition, constrained
        ** to solutions lower than 0.5
        */
        CN_leaf = 1.0 / s->shootnc;
        arg1 = MAX(0.0, (CN_leaf - CN_ref) / CN_ref);
        
        CP_leaf = 1.0 / s->shootpc;
        arg2 = MAX(0.0, (CP_leaf - CP_ref) / CP_ref);
        
        frac_to_rexc = MIN(0.5, MIN(p->a0rhizo + p->a1rhizo * arg1, p->a0rhizo + p->a1rhizo * arg2));
    }

    /* Rhizodeposition */
    f->root_exc = frac_to_rexc * f->cproot;
    if (float_eq(f->cproot, 0.0)) {
        f->root_exn = 0.0;
        f->root_exp = 0.0;
    } else {
        /*
        ** N flux associated with rhizodeposition is based on the assumption
        ** that the CN ratio of rhizodeposition is equal to that of fine root
        ** growth
        */
        f->root_exn = f->root_exc * (f->nproot / f->cproot);
        f->root_exp = f->root_exc * (f->pproot / f->cproot);
    }

    /*
    ** Need to remove exudation C & N fluxes from fine root growth fluxes so
    ** that things balance.
    */
    f->cproot -= f->root_exc;
    f->nproot -= f->root_exn;
    f->pproot -= f->root_exp;

    return;
}

void carbon_daily_production(control *c, fluxes *f, met *m, params *p, state *s,
                             double daylen) {
    /* Calculate GPP, NPP and plant respiration at the daily timestep

    Parameters:
    -----------
    daylen : float
        daytime length (hrs)

    References:
    -----------
    * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.
    */
    double leafn, leafp, fc, ncontent, pcontent;

    if (s->lai > 0.0) {
        /* average leaf nitrogen content (g N m-2 leaf) */
        leafn = (s->shootnc * p->cfracts / p->sla * KG_AS_G);
        /* average leaf phosphorus content (g P m-2 leaf) */
        leafp = (s->shootpc * p->cfracts / p->sla * KG_AS_G);
      
        /* total nitrogen content of the canopy */
        ncontent = leafn * s->lai;
        /* total phosphorus content of the canopy */
        pcontent = leafp * s->lai;

    } else {
        ncontent = 0.0;
        pcontent = 0.0;
    }

    /* When canopy is not closed, canopy light interception is reduced
        - calculate the fractional ground cover */
    if (s->lai < p->lai_closed) {
        /* discontinuous canopies */
        fc = s->lai / p->lai_closed;
    } else {
        fc = 1.0;
    }

    /* fIPAR - the fraction of intercepted PAR = IPAR/PAR incident at the
       top of the canopy, accounting for partial closure based on Jackson
       and Palmer (1979). */
    if (s->lai > 0.0)
        s->fipar = ((1.0 - exp(-p->kext * s->lai / fc)) * fc);
    else
        s->fipar = 0.0;

    if (c->water_stress) {
        /* Calculate the soil moisture availability factors [0,1] in the
           topsoil and the entire root zone */
        calculate_soil_water_fac(c, p, s);
    } else {
        /* really this should only be a debugging option! */
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }
    /* Estimate photosynthesis */
    if (c->assim_model == BEWDY){
        exit(EXIT_FAILURE);
    } else if (c->assim_model == MATE) {
        if (c->ps_pathway == C3) {
            mate_C3_photosynthesis(c, f, m, p, s, daylen, ncontent, pcontent);
        } else {
            mate_C4_photosynthesis(c, f, m, p, s, daylen, ncontent, pcontent);
        }
    } else {
        fprintf(stderr,"Unknown photosynthesis model'");
        exit(EXIT_FAILURE);
    }


    /* Calculate plant respiration */
    if (c->respiration_model == FIXED) {
        /* Plant respiration assuming carbon-use efficiency. */
        f->auto_resp = f->gpp * p->cue;
    } else if(c->respiration_model == TEMPERATURE) {
        fprintf(stderr, "Not implemented yet");
        exit(EXIT_FAILURE);
    } else if (c->respiration_model == BIOMASS) {
        fprintf(stderr, "Not implemented yet");
        exit(EXIT_FAILURE);
    }

    /* Calculate NPP */
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;

    return;
}

void calculate_cnp_wood_ratios(control *c, params *p, state *s, double npitfac, 
                             double nitfac, double pitfac, 
                             double *ncbnew, double *nccnew, 
                             double *ncwimm, double *ncwnew,
                             double *pcbnew, double *pccnew,
                             double *pcwimm, double *pcwnew) {
    /* Estimate the N:C and P:C ratio in the branch and stem. Option to vary
    the N:C and P:C ratio of the stem following Jeffreys (1999) or keep it a fixed
    fraction

    Parameters:
    -----------
    npitfac : float
        min of nitfac and pitfac;
    nitfac: float
       leaf N:C as a fraction of Ncmaxyoung;
    pitfac : float
       leaf P:C as a fraction of Pcmaxyoung;

    Returns:
    --------
    ncbnew : float
        N:C ratio of branch
    nccnew : double
        N:C ratio of coarse root
    ncwimm : float
        N:C ratio of immobile stem
    ncwnew : float
        N:C ratio of mobile stem
    pcbnew : float
        P:C ratio of branch
    pccnew : double
        P:C ratio of coarse root
    pcwimm : float
        P:C ratio of immobile stem
    pcwnew : float
        P:C ratio of mobile stem
     
    References:
    ----------
    * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
      with modelled implications for forest productivity under elevated
      atmospheric carbon dioxide. PhD.
    */

    /* calculate N:C ratios */
    if (npitfac < nitfac) {
      /* n:c ratio of new branch wood*/
      *ncbnew = p->ncbnew + pitfac * (p->ncbnew - p->ncbnewz);

      /* n:c ratio of coarse root */
      *nccnew = p->nccnew + pitfac * (p->nccnew - p->nccnewz);
      
      /* fixed N:C in the stemwood */
      if (c->fixed_stem_nc) {
        /* n:c ratio of stemwood - immobile pool and new ring */
        *ncwimm = p->ncwimm + pitfac * (p->ncwimm - p->ncwimmz);
        
        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = p->ncwnew + pitfac * (p->ncwnew - p->ncwnewz);
        
        /* vary stem N:C based on reln with foliage, see Jeffreys PhD thesis.
        Jeffreys 1999 showed that N:C ratio of new wood increases with foliar N:C ratio,
        modelled here based on evidence as a linear function. */
      } else {
        *ncwimm = MAX(0.0, (0.0282 * s->shootnc + 0.000234) * p->fhw);
        
        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = MAX(0.0, 0.162 * s->shootnc - 0.00143);
      }
    } else {
      /* n:c ratio of new branch wood*/
      *ncbnew = p->ncbnew + nitfac * (p->ncbnew - p->ncbnewz);
      
      /* n:c ratio of coarse root */
      *nccnew = p->nccnew + nitfac * (p->nccnew - p->nccnewz);
      
      /* fixed N:C in the stemwood */
      if (c->fixed_stem_nc) {
        /* n:c ratio of stemwood - immobile pool and new ring */
        *ncwimm = p->ncwimm + nitfac * (p->ncwimm - p->ncwimmz);
        
        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = p->ncwnew + nitfac * (p->ncwnew - p->ncwnewz);
        
        /* vary stem N:C based on reln with foliage, see Jeffreys. Jeffreys 1999
        showed that N:C ratio of new wood increases with foliar N:C ratio,
        modelled here based on evidence as a linear function. */
      } else {
        *ncwimm = MAX(0.0, (0.0282 * s->shootnc + 0.000234) * p->fhw);
        
        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = MAX(0.0, 0.162 * s->shootnc - 0.00143);
      }
    }
    
    /* calculate P:C ratios */
    if (npitfac < nitfac) {
      /* p:c ratio of new branch wood*/
      *pcbnew = p->pcbnew + pitfac * (p->pcbnew - p->pcbnewz);
      
      /* p:c ratio of coarse root */
      *pccnew = p->pccnew + pitfac * (p->pccnew - p->pccnewz);
      
      /* fixed P:C in the stemwood */
      if (c->fixed_stem_pc) {
        /* p:c ratio of stemwood - immobile pool and new ring */
        *pcwimm = p->pcwimm + pitfac * (p->pcwimm - p->pcwimmz);
        
        /* New stem ring P:C at critical leaf P:C (mobile) */
        *pcwnew = p->pcwnew + pitfac * (p->pcwnew - p->pcwnewz);
        
        /* vary stem P:C based on reln with foliage, see Jeffreys. Jeffreys 1999
        showed that P:C ratio of new wood increases with foliar P:C ratio,
        modelled here based on evidence as a linear function. */
      } else {
        *pcwimm = MAX(0.0, (0.0282 * s->shootpc + 0.000234) * p->fhwp);
        
        /* New stem ring P:C at critical leaf P:C (mobile) */
        *pcwnew = MAX(0.0, 0.162 * s->shootpc - 0.00143);
      }
    } else {
      /* p:c ratio of new branch wood*/
      *pcbnew = p->pcbnew + nitfac * (p->pcbnew - p->pcbnewz);
      
      /* p:c ratio of coarse root */
      *pccnew = p->pccnew + nitfac * (p->pccnew - p->pccnewz);
      
      /* fixed P:C in the stemwood */
      if (c->fixed_stem_pc) {
        /* p:c ratio of stemwood - immobile pool and new ring */
        *pcwimm = p->pcwimm + nitfac * (p->pcwimm - p->pcwimmz);
        
        /* New stem ring P:C at critical leaf P:C (mobile) */
        *pcwnew = p->pcwnew + nitfac * (p->pcwnew - p->pcwnewz);
        
        /* vary stem P:C based on reln with foliage, see Jeffreys. Jeffreys 1999
        showed that P:C ratio of new wood increases with foliar P:C ratio,
        modelled here based on evidence as a linear function. */
      } else {
        *pcwimm = MAX(0.0, (0.0282 * s->shootpc + 0.000234) * p->fhwp);
        
        /* New stem ring P:C at critical leaf P:C (mobile) */
        *pcwnew = MAX(0.0, 0.162 * s->shootpc - 0.00143);
      }
    }
    
    
    return;
}



int np_allocation(control *c, fluxes *f, params *p, state *s,
                        double ncbnew, double nccnew, double ncwimm,
                        double ncwnew, 
                        double pcbnew, double pccnew, double pcwimm,
                        double pcwnew,double fdecay, double rdecay, int doy) {
    /* Nitrogen and phosphorus distribution - allocate available N and P (mineral) through system.
    N and P is first allocated to the woody component, surplus N and P is then allocated
    to the shoot and roots with flexible ratios.

    References:
    -----------
    McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.

    Parameters:
    -----------
    ncbnew : float
        N:C ratio of branch
    ncwimm : float
        N:C ratio of immobile stem
    ncwnew : float
        N:C ratio of mobile stem
    pcbnew : float
        P:C ratio of branch
    pcwimm : float
        P:C ratio of immobile stem
    pcwnew : float
        P:C ratio of mobile stem
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate
    */

    int    recalc_wb;
    double nsupply, psupply, rtot, ntot, ptot, arg1, arg2, lai_inc = 0.0, conv;
    double depth_guess = 1.0;

    /* default is we don't need to recalculate the water balance,
       however if we cut back on NPP due to available N and P below then we do
       need to do this */
    recalc_wb = FALSE;

    /* N and P retranslocated proportion from dying plant tissue and stored within
       the plant */
    f->retrans = nitrogen_retrans(c, f, p, s, fdecay, rdecay, doy);
    f->retransp = phosphorus_retrans(c, f, p, s, fdecay, rdecay, doy);
    f->nuptake = calculate_nuptake(c, p, s);
    f->puptake = calculate_puptake(c, p, s);
    
    /* diagnosis */ 
    
    fprintf(stderr, "puptake = %f\n", f->puptake);
    
     
    /*  Ross's Root Model. */
    if (c->model_optroot) {

        /* convert t ha-1 day-1 to gN m-2 year-1 */
        nsupply = (calculate_nuptake(c, p, s) *
                   TONNES_HA_2_G_M2 * DAYS_IN_YRS);
      
        psupply = (calculate_puptake(c, p, s) *
                  TONNES_HA_2_G_M2 * DAYS_IN_YRS);

        /* covnert t ha-1 to kg DM m-2 */
        rtot = s->root * TONNES_HA_2_KG_M2 / p->cfracts;
        /*f->nuptake_old = f->nuptake; */

        calc_opt_root_depth(p->d0x, p->r0, p->topsoil_depth * MM_TO_M,
                            rtot, nsupply, psupply, depth_guess, &s->root_depth,
                            &f->nuptake, &f->puptake, &f->rabove);

        /*umax = self.rm.calc_umax(f->nuptake) */

        /* covert nuptake from gN m-2 year-1  to t ha-1 day-1 */
        f->nuptake = f->nuptake * G_M2_2_TONNES_HA * YRS_IN_DAYS;
        f->puptake = f->puptake * G_M2_2_TONNES_HA * YRS_IN_DAYS;
        

        /* covert from kg DM N m-2 to t ha-1 */
        f->deadroots = p->rdecay * f->rabove * p->cfracts * KG_M2_2_TONNES_HA;
        f->deadrootn = s->rootnc * (1.0 - p->rretrans) * f->deadroots;
        f->deadrootp = s->rootpc * (1.0 - p->rretrans) * f->deadroots;
        
    }

    /* Mineralised nitrogen lost from the system by volatilisation/leaching */
    f->nloss = p->rateloss * s->inorgn;
    
    /* Mineralised P lost from the system by leaching */
    f->ploss = p->prateloss * s->inorglabp;

    /* total nitrogen/phosphorus to allocate */
    ntot = MAX(0.0, f->nuptake + f->retrans);
    ptot = MAX(0.0, f->puptake + f->retransp);
    

    if (c->deciduous_model) {
        /* allocate N to pools with fixed N:C ratios */

        /* N flux into new ring (immobile component -> structrual components) */
        f->npstemimm = f->wnimrate * s->growing_days[doy];

        /* N flux into new ring (mobile component -> can be retrans for new
           woody tissue) */
        f->npstemmob = f->wnmobrate * s->growing_days[doy];
        f->nproot = s->n_to_alloc_root / c->num_days;
        f->npcroot = f->cnrate * s->growing_days[doy];
        f->npleaf = f->lnrate * s->growing_days[doy];
        f->npbranch = f->bnrate * s->growing_days[doy];
        
        /* allocate P to pools with fixed P:C ratios */
        f->ppstemimm = f->wpimrate * s->growing_days[doy];
        f->ppstemmob = f->wpmobrate * s->growing_days[doy];
        f->pproot = s->p_to_alloc_root / c->num_days;
        f->ppcroot = f->cprate * s->growing_days[doy];
        f->ppleaf = f->lprate * s->growing_days[doy];
        f->ppbranch = f->bprate * s->growing_days[doy];
        
    } else {
        /* allocate N to pools with fixed N:C ratios */

        /* N flux into new ring (immobile component -> structural components) */
        f->npstemimm = f->npp * f->alstem * ncwimm;

        /* N flux into new ring (mobile component -> can be retrans for new
           woody tissue) */
        f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
        f->npbranch = f->npp * f->albranch * ncbnew;
        f->npcroot = f->npp * f->alcroot * nccnew;

        /* allocate P to pools with fixed P:C ratios */
        f->ppstemimm = f->npp * f->alstem * pcwimm;
        f->ppstemmob = f->npp * f->alstem * (pcwnew - pcwimm);
        f->ppbranch = f->npp * f->albranch * pcbnew;
        f->ppcroot = f->npp * f->alcroot * pccnew;
        
        
        /* If we have allocated more N than we have available
            - cut back C prodn */
        arg1 = f->npstemimm + f->npstemmob + f->npbranch + f->npcroot;

        /* If we have allocated more P than we have available
         - cut back C prodn */
        arg2 = f->ppstemimm + f->ppstemmob + f->ppbranch + f->ppcroot;
        
        /* diagnosis */ 
        /*
        fprintf(stderr, "arg1 = %f\n", arg1);
        fprintf(stderr, "ntot = %f\n", ntot);
        fprintf(stderr, "arg2 = %f\n", arg2);
        fprintf(stderr, "ptot = %f\n", ptot);
        */

        if (arg1 > ntot && c->fixleafnc == FALSE && c->fixed_lai && c->ncycle) {

            /* Need to readjust the LAI for the reduced growth as this will
               have already been increased. First we need to figure out how
               much we have increased LAI by, important it is done here
               before cpleaf is reduced! */
            if (float_eq(s->shoot, 0.0)) {
                lai_inc = 0.0;
            } else {
                lai_inc = (f->cpleaf *
                           (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
                           (f->deadleaves + f->ceaten) * s->lai / s->shoot);
            }

            f->npp *= ntot / (f->npstemimm + f->npstemmob + \
                              f->npbranch + f->npcroot);

            /* need to adjust growth values accordingly as well */
            f->cpleaf = f->npp * f->alleaf;
            f->cproot = f->npp * f->alroot;
            f->cpcroot = f->npp * f->alcroot;
            f->cpbranch = f->npp * f->albranch;
            f->cpstem = f->npp * f->alstem;

            f->npbranch = f->npp * f->albranch * ncbnew;
            f->npstemimm = f->npp * f->alstem * ncwimm;
            f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
            f->npcroot = f->npp * f->alcroot * nccnew;
            
            f->ppbranch = f->npp * f->albranch * pcbnew;
            f->ppstemimm = f->npp * f->alstem * pcwimm;
            f->ppstemmob = f->npp * f->alstem * (pcwnew - pcwimm);
            f->ppcroot = f->npp * f->alcroot * pccnew;

            /* Save WUE before cut back */
            f->wue = f->gpp_gCm2 / f->transpiration;

            /* Also need to recalculate GPP and thus Ra and return a flag
               so that we know to recalculate the water balance. */
            f->gpp = f->npp / p->cue;
            conv = G_AS_TONNES / M2_AS_HA;
            f->gpp_gCm2 = f->gpp / conv;
            f->gpp_am = f->gpp_gCm2 / 2.0;
            f->gpp_pm = f->gpp_gCm2 / 2.0;

            /* New respiration flux */
            f->auto_resp =  f->gpp - f->npp;
            recalc_wb = TRUE;

            /* Now reduce LAI for down-regulated growth. */
            if (c->deciduous_model) {
                if (float_eq(s->shoot, 0.0)) {
                    s->lai = 0.0;
                } else if (s->leaf_out_days[doy] > 0.0) {
                    s->lai -= lai_inc;
                    s->lai += (f->cpleaf *
                               (p->sla * M2_AS_HA / \
                               (KG_AS_TONNES * p->cfracts)) -
                               (f->deadleaves + f->ceaten) * s->lai / s->shoot);
                } else {
                    s->lai = 0.0;
                }
            } else {
                /* update leaf area [m2 m-2] */
                if (float_eq(s->shoot, 0.0)) {
                    s->lai = 0.0;
                } else {
                    s->lai -= lai_inc;
                    s->lai += (f->cpleaf *
                               (p->sla * M2_AS_HA / \
                               (KG_AS_TONNES * p->cfracts)) -
                               (f->deadleaves + f->ceaten) * s->lai / s->shoot);
                }
            }
        }
        
        if (arg2 > ptot && c->fixleafpc == FALSE && c->fixed_lai && c->pcycle) {
          
          /* Need to readjust the LAI for the reduced growth as this will
          have already been increased. First we need to figure out how
          much we have increased LAI by, important it is done here
          before cpleaf is reduced! */
          if (float_eq(s->shoot, 0.0)) {
            lai_inc = 0.0;
          } else {
            lai_inc = (f->cpleaf *
              (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
              (f->deadleaves + f->ceaten) * s->lai / s->shoot);
          }
          
          f->npp *= ptot / (f->ppstemimm + f->ppstemmob + \
            f->ppbranch + f->ppcroot);
          
          /* need to adjust growth values accordingly as well */
          f->cpleaf = f->npp * f->alleaf;
          f->cproot = f->npp * f->alroot;
          f->cpcroot = f->npp * f->alcroot;
          f->cpbranch = f->npp * f->albranch;
          f->cpstem = f->npp * f->alstem;
          
          f->npbranch = f->npp * f->albranch * ncbnew;
          f->npstemimm = f->npp * f->alstem * ncwimm;
          f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
          f->npcroot = f->npp * f->alcroot * nccnew;
          
          f->ppbranch = f->npp * f->albranch * pcbnew;
          f->ppstemimm = f->npp * f->alstem * pcwimm;
          f->ppstemmob = f->npp * f->alstem * (pcwnew - pcwimm);
          f->ppcroot = f->npp * f->alcroot * pccnew;
          
          /* Save WUE before cut back */
          f->wue = f->gpp_gCm2 / f->transpiration;
          
          /* Also need to recalculate GPP and thus Ra and return a flag
          so that we know to recalculate the water balance. */
          f->gpp = f->npp / p->cue;
          conv = G_AS_TONNES / M2_AS_HA;
          f->gpp_gCm2 = f->gpp / conv;
          f->gpp_am = f->gpp_gCm2 / 2.0;
          f->gpp_pm = f->gpp_gCm2 / 2.0;
          
          
          /* New respiration flux */
          f->auto_resp =  f->gpp - f->npp;
          recalc_wb = TRUE;
          
          /* Now reduce LAI for down-regulated growth. */
          if (c->deciduous_model) {
            if (float_eq(s->shoot, 0.0)) {
              s->lai = 0.0;
            } else if (s->leaf_out_days[doy] > 0.0) {
              s->lai -= lai_inc;
              s->lai += (f->cpleaf *
                (p->sla * M2_AS_HA / \
                (KG_AS_TONNES * p->cfracts)) -
                (f->deadleaves + f->ceaten) * s->lai / s->shoot);
            } else {
              s->lai = 0.0;
            }
          } else {
            /* update leaf area [m2 m-2] */
            if (float_eq(s->shoot, 0.0)) {
              s->lai = 0.0;
            } else {
              s->lai -= lai_inc;
              s->lai += (f->cpleaf *
                (p->sla * M2_AS_HA / \
                (KG_AS_TONNES * p->cfracts)) -
                (f->deadleaves + f->ceaten) * s->lai / s->shoot);
            }
          }
        }

        /* Nitrogen reallocation to flexible-ratio pools */
        ntot -= f->npbranch + f->npstemimm + f->npstemmob + f->npcroot;
        ntot = MAX(0.0, ntot);

        /* allocate remaining N to flexible-ratio pools */
        f->npleaf = ntot * f->alleaf / (f->alleaf + f->alroot * p->ncrfac);
        f->nproot = ntot - f->npleaf;
        
        /* Phosphorus reallocation to flexible-ratio pools */
        ptot -= f->ppbranch + f->ppstemimm + f->ppstemmob + f->ppcroot;
        ptot = MAX(0.0, ptot);
        
        /* allocate remaining P to flexible-ratio pools */
        f->ppleaf = ptot * f->alleaf / (f->alleaf + f->alroot * p->pcrfac);
        f->pproot = ptot - f->ppleaf;
        
    }
    return (recalc_wb);
}


double calculate_growth_stress_limitation(params *p, state *s) {
    /* Calculate level of stress due to nitrogen, phosphorus or water availability */
    double nlim, plim, current_limitation;
    double nc_opt = 0.04;
    double pc_opt = 0.004;

    /* N limitation based on leaf NC ratio */
    if (s->shootnc < p->nf_min) {
        nlim = 0.0;
    } else if (s->shootnc < nc_opt && s->shootnc > p->nf_min) {
        nlim = 1.0 - ((nc_opt - s->shootnc) / (nc_opt - p->nf_min));
    } else {
        nlim = 1.0;
    }
    
    /* P limitation based on leaf PC ratio */
    if (s->shootpc < p->pf_min) {
      plim = 0.0;
    } else if (s->shootpc < pc_opt && s->shootpc > p->pf_min) {
      plim = 1.0 - ((pc_opt - s->shootpc) / (pc_opt - p->pf_min));
    } else {
      plim = 1.0;
    }

    /*
     * Limitation by nitrogen, water and phosphorus. Water constraint is implicit,
     * in that, water stress results in an increase of root mass,
     * which are assumed to spread horizontally within the rooting zone.
     * So in effect, building additional root mass doesnt alleviate the
     * water limitation within the model. However, it does more
     * accurately reflect an increase in root C production at a water
     * limited site. This implementation is also consistent with other
     * approaches, e.g. LPJ. In fact I dont see much evidence for models
     * that have a flexible bucket depth. Minimum constraint is limited to
     * 0.1, following Zaehle et al. 2010 (supp), eqn 18.
     */
    current_limitation = MAX(0.1, MIN(nlim, MIN(plim, s->wtfac_root)));
    return (current_limitation);
}


void calc_carbon_allocation_fracs(control *c, fluxes *f, params *p, state *s,
                                  double npitfac) {
    /* Carbon allocation fractions to move photosynthate through the plant.

    Parameters:
    -----------
    npitfac : float
        the smallest value of leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0) &
        leaf P:C as a fraction of "Pcmaxfyoung" (max 1.0)

    Returns:
    --------
    alleaf : float
        allocation fraction for shoot
    alroot : float
        allocation fraction for fine roots
    albranch : float
        allocation fraction for branches
    alstem : float
        allocation fraction for stem

    References:
    -----------
    Corbeels, M. et al (2005) Ecological Modelling, 187, 449-474.
    McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
    */
    double min_leaf_alloc, adj, arg1, arg2, arg3, arg4, leaf2sa_target,
           sap_cross_sec_area, lr_max, stress, mis_match, orig_ar,
           reduction, target_branch, coarse_root_target, left_over,
           total_alloc, leaf2sap, spare;

    /* this is obviously arbitary */
    double min_stem_alloc = 0.01;

    if (c->alloc_model == FIXED){
        f->alleaf = (p->c_alloc_fmax + npitfac *
                     (p->c_alloc_fmax - p->c_alloc_fmin));

        f->alroot = (p->c_alloc_rmax + npitfac *
                     (p->c_alloc_rmax - p->c_alloc_rmin));

        f->albranch = (p->c_alloc_bmax + npitfac *
                       (p->c_alloc_bmax - p->c_alloc_bmin));

        /* allocate remainder to stem */
        f->alstem = 1.0 - f->alleaf - f->alroot - f->albranch;

        f->alcroot = p->c_alloc_cmax * f->alstem;
        f->alstem -= f->alcroot;

    } else if (c->alloc_model == GRASSES) {

        /* First figure out root allocation given available water & nutrients
           hyperbola shape to allocation */
        f->alroot = (p->c_alloc_rmax * p->c_alloc_rmin /
                     (p->c_alloc_rmin + (p->c_alloc_rmax - p->c_alloc_rmin) *
                      s->prev_sma));
        f->alleaf = 1.0 - f->alroot;

        /* Now adjust root & leaf allocation to maintain balance, accounting
           for stress e.g. -> Sitch et al. 2003, GCB.

         leaf-to-root ratio under non-stressed conditons
        lr_max = 0.8;

         Calculate adjustment on lr_max, based on current "stress"
           calculated from running mean of N and water stress
        stress = lr_max * s->prev_sma;

        calculate new allocation fractions based on imbalance in *biomass*
        mis_match = s->shoot / (s->root * stress);


        if (mis_match > 1.0) {
            reduce leaf allocation fraction
            adj = f->alleaf / mis_match;
            f->alleaf = MAX(p->c_alloc_fmin, MIN(p->c_alloc_fmax, adj));
            f->alroot = 1.0 - f->alleaf;
        } else {
             reduce root allocation
            adj = f->alroot * mis_match;
            f->alroot = MAX(p->c_alloc_rmin, MIN(p->c_alloc_rmax, adj));
            f->alleaf = 1.0 - f->alroot;
        }*/
        f->alstem = 0.0;
        f->albranch = 0.0;
        f->alcroot = 0.0;

    } else if (c->alloc_model == ALLOMETRIC) {

        /* Calculate tree height: allometric reln using the power function
           (Causton, 1985) */
        s->canht = p->heighto * pow(s->stem, p->htpower);

        /* LAI to stem sapwood cross-sectional area (As m-2 m-2)
           (dimensionless)
           Assume it varies between LS0 and LS1 as a linear function of tree
           height (m) */
        arg1 = s->sapwood * TONNES_AS_KG * M2_AS_HA;
        arg2 = s->canht * p->density * p->cfracts;
        sap_cross_sec_area = arg1 / arg2;
        leaf2sap = s->lai / sap_cross_sec_area;

        /* Allocation to leaves dependant on height. Modification of pipe
           theory, leaf-to-sapwood ratio is not constant above a certain
           height, due to hydraulic constraints (Magnani et al 2000; Deckmyn
           et al. 2006). */

        if (s->canht < p->height0) {
            leaf2sa_target = p->leafsap0;
        } else if (float_eq(s->canht, p->height1)) {
            leaf2sa_target = p->leafsap1;
        } else if (s->canht > p->height1) {
            leaf2sa_target = p->leafsap1;
        } else {
            arg1 = p->leafsap0;
            arg2 = p->leafsap1 - p->leafsap0;
            arg3 = s->canht - p->height0;
            arg4 = p->height1 - p->height0;
            leaf2sa_target = arg1 + (arg2 * arg3 / arg4);
        }
        f->alleaf = alloc_goal_seek(leaf2sap, leaf2sa_target, p->c_alloc_fmax,
                                    p->targ_sens);

        /* Allocation to branch dependent on relationship between the stem
           and branch */
        target_branch = p->branch0 * pow(s->stem, p->branch1);
        f->albranch = alloc_goal_seek(s->branch, target_branch, p->c_alloc_bmax,
                                      p->targ_sens);

        coarse_root_target = p->croot0 * pow(s->stem, p->croot1);
        f->alcroot = alloc_goal_seek(s->croot, coarse_root_target,
                                      p->c_alloc_cmax, p->targ_sens);

        /* figure out root allocation given available water & nutrients
           hyperbola shape to allocation, this is adjusted below as we aim
           to maintain a functional balance */

        f->alroot = (p->c_alloc_rmax * p->c_alloc_rmin /
                     (p->c_alloc_rmin + (p->c_alloc_rmax - p->c_alloc_rmin) *
                      s->prev_sma));

        f->alstem = 1.0 - f->alroot - f->albranch - f->alleaf - f->alcroot;

        /* minimum allocation to leaves - without it tree would die, as this
           is done annually. */
        if (c->deciduous_model) {
            if (f->alleaf < 0.05) {
                min_leaf_alloc = 0.05;
                if (f->alstem > min_leaf_alloc)
                    f->alstem -= min_leaf_alloc;
                else
                    f->alroot -= min_leaf_alloc;
                f->alleaf = min_leaf_alloc;
            }
        }
    } else {
        fprintf(stderr, "Unknown C allocation model: %d\n", c->alloc_model);
        exit(EXIT_FAILURE);
    }

    /*printf("%f %f %f %f %f\n", f->alleaf, f->albranch + f->alstem, f->alroot,  f->alcroot, s->canht);*/

    /* Total allocation should be one, if not print warning */
    total_alloc = f->alroot + f->alleaf + f->albranch + f->alstem + f->alcroot;
    if (total_alloc > 1.0+EPSILON) {
        fprintf(stderr, "Allocation fracs > 1: %.13f\n", total_alloc);
        exit(EXIT_FAILURE);
    }

    return;
}


double alloc_goal_seek(double simulated, double target, double alloc_max,
                       double sensitivity) {

    /* Sensitivity parameter characterises how allocation fraction respond
       when the leaf:sapwood area ratio departs from the target value
       If sensitivity close to 0 then the simulated leaf:sapwood area ratio
       will closely track the target value */
    double frac = 0.5 + 0.5 * (1.0 - simulated / target) / sensitivity;

    return MAX(0.0, alloc_max * MIN(1.0, frac));
}

void carbon_allocation(control *c, fluxes *f, params *p, state *s,
                       double npitfac, int doy) {
    /* C distribution - allocate available C through system

    Parameters:
    -----------
    npitfac : float
        leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0) and leaf P:C as a fraction of "Mcmaxfyoung"
    */
    double days_left;
    if (c->deciduous_model) {
        days_left = s->growing_days[doy];
        f->cpleaf = f->lrate * days_left;
        f->cpbranch = f->brate * days_left;
        f->cpstem = f->wrate * days_left;
        f->cproot = s->c_to_alloc_root * 1.0 / c->num_days;
        f->cpcroot = f->crate * days_left;
    } else {
        f->cpleaf = f->npp * f->alleaf;
        f->cproot = f->npp * f->alroot;
        f->cpcroot = f->npp * f->alcroot;
        f->cpbranch = f->npp * f->albranch;
        f->cpstem = f->npp * f->alstem;
    }

    /* evaluate SLA of new foliage accounting for variation in SLA
       with tree and leaf age (Sands and Landsberg, 2002). Assume
       SLA of new foliage is linearly related to leaf N:C ratio
       via nitfac. Based on date from two E.globulus stands in SW Aus, see
       Corbeels et al (2005) Ecological Modelling, 187, 449-474.
       (m2 onesided/kg DW) */
    p->sla = p->slazero + npitfac * (p->slamax - p->slazero);

    if (c->deciduous_model) {
        if (float_eq(s->shoot, 0.0)) {
            s->lai = 0.0;
        } else if (s->leaf_out_days[doy] > 0.0) {
            s->lai += (f->cpleaf *
                      (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
                      (f->deadleaves + f->ceaten) * s->lai / s->shoot);
        } else {
            s->lai = 0.0;
        }
    } else {
        /* update leaf area [m2 m-2] */
        if (float_eq(s->shoot, 0.0)) {
            s->lai = 0.0;
        } else {
            s->lai += (f->cpleaf *
                      (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
                      (f->deadleaves + f->ceaten) * s->lai / s->shoot);
        }
    }

    if (c->fixed_lai) {
        s->lai = p->fix_lai;
    }

    return;
}

void update_plant_state(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
    /*
    Daily change in C content

    Parameters:
    -----------
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate

    */

    double age_effect;
    double ncmaxf, ncmaxr;
    double extrasn, extrarn;   /* extra_s_n and extra_r_n - extra shoot/root n uptake */
    double pcmaxf, pcmaxr;
    double extrasp, extrarp;   /* extra_s_p and extra_r_p - extra shoot/root p uptake */

    /*
    ** Carbon pools
    */
    s->shoot += f->cpleaf - f->deadleaves - f->ceaten;
    s->root += f->cproot - f->deadroots;
    s->croot += f->cpcroot - f->deadcroots;
    s->branch += f->cpbranch - f->deadbranch;
    s->stem += f->cpstem - f->deadstems;

    /* annoying but can't see an easier way with the code as it is.
       If we are modelling grases, i.e. no stem them without this
       the sapwood will end up being reduced to a silly number as
       deadsapwood will keep being removed from the pool, even though there
       is no wood. */
    if (float_eq(s->stem, 0.01)) {
        s->sapwood = 0.01;
    } else if (s->stem < 0.01) {
        s->sapwood = 0.01;
    } else {
        s->sapwood += f->cpstem - f->deadsapwood;
    }


    /*
    ** Nitrogen and Phosphorus pools
    */
    if (c->deciduous_model) {
        s->shootn += (f->npleaf - (f->lnrate * s->remaining_days[doy]) -
                      f->neaten);
        s->shootp += (f->ppleaf - (f->lprate * s->remaining_days[doy]) -
                      f->peaten);
    } else {
        s->shootn += f->npleaf - fdecay * s->shootn - f->neaten;
        s->shootp += f->ppleaf - fdecay * s->shootp - f->peaten;
    }

    s->branchn += f->npbranch - p->bdecay * s->branchn;
    s->rootn += f->nproot - rdecay * s->rootn;
    s->crootn += f->npcroot - p->crdecay * s->crootn;
    s->stemnimm += f->npstemimm - p->wdecay * s->stemnimm;
    s->stemnmob += (f->npstemmob - p->wdecay * s->stemnmob - p->retransmob *
                    s->stemnmob);
    s->stemn = s->stemnimm + s->stemnmob;
    
    s->branchp += f->ppbranch - p->bdecay * s->branchp;
    s->rootp += f->pproot - rdecay * s->rootp;
    s->crootp += f->ppcroot - p->crdecay * s->crootp;
    s->stempimm += f->ppstemimm - p->wdecay * s->stempimm;
    s->stempmob += (f->ppstemmob - p->wdecay * s->stempmob - p->retransmob *
                    s->stempmob);
    s->stemp = s->stempimm + s->stempmob;


    if (c->deciduous_model == FALSE) {
        /*
           =============================
            Enforce maximum N:C and P:C ratios.
           =============================
        */

        /* If foliage or root N/C exceeds its max, then N uptake is cut back
           Similarly, of foliage or root P/C exceeds max, then P uptake is cut back */

        /* maximum leaf n:c and p:c ratios is function of stand age
            - switch off age effect by setting ncmaxfyoung = ncmaxfold 
            - switch off age effect by setting pcmaxfyoung = pcmaxfold*/
        age_effect = (s->age - p->ageyoung) / (p->ageold - p->ageyoung);
        ncmaxf = p->ncmaxfyoung - (p->ncmaxfyoung - p->ncmaxfold) * age_effect;
        pcmaxf = p->pcmaxfyoung - (p->pcmaxfyoung - p->pcmaxfold) * age_effect;

        if (ncmaxf < p->ncmaxfold)
            ncmaxf = p->ncmaxfold;

        if (ncmaxf > p->ncmaxfyoung)
            ncmaxf = p->ncmaxfyoung;
        
        if (pcmaxf < p->pcmaxfold)
            pcmaxf = p->pcmaxfold;
        
        if (pcmaxf > p->pcmaxfyoung)
            pcmaxf = p->pcmaxfyoung;

        extrasn = 0.0;
        if (s->lai > 0.0) {

            if (s->shootn > (s->shoot * ncmaxf)) {
                extrasn = s->shootn - s->shoot * ncmaxf;

                /* Ensure N uptake cannot be reduced below zero. */
                if (extrasn >  f->nuptake)
                    extrasn = f->nuptake;

                s->shootn -= extrasn;
                f->nuptake -= extrasn;
            }
        }
        
        extrasp = 0.0;
        if (s->lai > 0.0) {
          
          if (s->shootp > (s->shoot * pcmaxf)) {
            extrasp = s->shootp - s->shoot * pcmaxf;
            
            /* Ensure P uptake cannot be reduced below zero. */
            if (extrasp >  f->puptake)
              extrasp = f->puptake;
            
            s->shootp -= extrasp;
            f->puptake -= extrasp;
          }
        }

        /* if root N:C ratio exceeds its max, then nitrogen uptake is cut
           back. n.b. new ring n/c max is already set because it is related
           to leaf n:c */

        /* max root n:c */
        ncmaxr = ncmaxf * p->ncrfac;
        extrarn = 0.0;
        if (s->rootn > (s->root * ncmaxr)) {
            extrarn = s->rootn - s->root * ncmaxr;

            /* Ensure N uptake cannot be reduced below zero. */
            if ((extrasn + extrarn) > f->nuptake)
                extrarn = f->nuptake - extrasn;

            s->rootn -= extrarn;
            f->nuptake -= extrarn;
        }

        /* max root p:c */
        pcmaxr = pcmaxf * p->pcrfac;
        extrarp = 0.0;
        if (s->rootp > (s->root * pcmaxr)) {
          extrarp = s->rootp - s->root * pcmaxr;
          
          /* Ensure P uptake cannot be reduced below zero. */
          if ((extrasp + extrarp) > f->puptake)
            extrarp = f->puptake - extrasp;
          
          s->rootp -= extrarp;
          f->puptake -= extrarp;
        }
    }
    
    
    /* Update deciduous storage pools */
    if (c->deciduous_model)
        calculate_cnp_store(c, f, s);

    return;
}

void precision_control(fluxes *f, state *s) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-10;

    /* C, N & P state variables */
    if (s->shoot < tolerance) {
        f->deadleaves += s->shoot;
        f->deadleafn += s->shootn;
        f->deadleafp += s->shootp;
        s->shoot = 0.0;
        s->shootn = 0.0;   
        s->shootp = 0.0;
    }

    if (s->branch < tolerance) {
        f->deadbranch += s->branch;
        f->deadbranchn += s->branchn;
        f->deadbranchp += s->branchp;
        s->branch = 0.0;
        s->branchn = 0.0;
        s->branchp = 0.0;
    }

    if (s->root < tolerance) {
        f->deadrootn += s->rootn;
        f->deadrootp += s->rootp;
        f->deadroots += s->root;
        s->root = 0.0;
        s->rootn = 0.0;
        s->rootp = 0.0;
    }

    if (s->croot < tolerance) {
        f->deadcrootn += s->crootn;
        f->deadcrootp += s->crootp;
        f->deadcroots += s->croot;
        s->croot = 0.0;
        s->crootn = 0.0;
        s->crootp = 0.0;
    }

    /* Not setting these to zero as this just leads to errors with desert
       regrowth...instead seeding them to a small value with a CN~25 and CP~300. */

    if (s->stem < tolerance) {
        f->deadstems += s->stem;
        f->deadstemn += s->stemn;
        f->deadstemp += s->stemp;
        s->stem = 0.001;
        s->stemn = 0.00004;
        s->stemp = 0.000003;
        s->stemnimm = 0.00004;
        s->stemnmob = 0.0;
        s->stempimm = 0.000003;
        s->stempmob = 0.0;
    }

    /* need separate one as this will become very small if there is no
       mobile stem N/P */
    if (s->stemnmob < tolerance) {
        f->deadstemn += s->stemnmob;
        s->stemnmob = 0.0;
    }

    if (s->stemnimm < tolerance) {
        f->deadstemn += s->stemnimm;
        s->stemnimm = 0.00004;
    }
    
    if (s->stempmob < tolerance) {
      f->deadstemp += s->stempmob;
      s->stempmob = 0.0;
    }
    
    if (s->stempimm < tolerance) {
      f->deadstemp += s->stempimm;
      s->stempimm = 0.000003;
    }

    return;
}


void calculate_cnp_store(control *c, fluxes *f, state *s) {
    /*
    Calculate labile C, N and P stores from which growth is allocated in the
    following year.
    */


    s->cstore += f->npp;
    s->nstore += f->nuptake + f->retrans;
    s->pstore += f->puptake + f->retransp;
    s->anpp += f->npp;

    /*
    double nstore_max, excess, k;
    double leaf_nc_max = 0.04;
    double CN_max = 100.0;


    if (c->alloc_model == GRASSES) {
        k = 0.3;
        nstore_max = MAX(1E-04, k * s->root * leaf_nc_max);
    } else {
        k = 0.15;
        nstore_max = MAX(1E-04, k * s->sapwood * leaf_nc_max);
    }

    s->nstore += f->nuptake + f->retrans;
    if (s->nstore > nstore_max) {
        s->nstore = nstore_max;
        f->nuptake = 0.0;
        f->retrans = 0.0;
    }

    s->cstore += f->npp;
    */
    /*
    if (s->cstore/nstore_max > CN_max) {
        excess = s->cstore - (nstore_max * CN_max);
        s->cstore = nstore_max * CN_max;
        f->auto_resp += excess;
    }*/


    return;
}


void calculate_average_alloc_fractions(fluxes *f, state *s, int days) {
    double excess;

    s->avg_alleaf /= (float) days;
    s->avg_alroot /= (float) days;
    s->avg_alcroot /= (float) days;
    s->avg_albranch /= (float) days;
    s->avg_alstem /= (float) days;

    f->alleaf = s->avg_alleaf;
    f->alroot = s->avg_alroot;
    f->alcroot = s->avg_alcroot;
    f->albranch = s->avg_albranch;
    f->alstem = s->avg_alstem;

    /*
        Because we are taking the average growing season fracs the total may
        end up being just under 1, due to rounding. So put the missing
        bit into the leaves - arbitary decision there
    */
    excess = 1.0 - f->alleaf - f->alroot - f->alcroot - f->albranch - f->alstem;
    f->alleaf += excess;

    return;
}

void allocate_stored_cnp(fluxes *f, params *p, state *s) {
    /*
    Allocate stored C,N and P. This is either down as the model is initialised
    for the first time or at the end of each year.
    */
    double ntot, ptot;

    /* ========================
       Carbon - fixed fractions
       ======================== */
    s->c_to_alloc_shoot = f->alleaf * s->cstore;
    s->c_to_alloc_root = f->alroot * s->cstore;
    s->c_to_alloc_croot = f->alcroot * s->cstore;
    s->c_to_alloc_branch = f->albranch * s->cstore;
    s->c_to_alloc_stem = f->alstem * s->cstore;

    /* =========================================================
        Nitrogen - Fixed ratios N allocation to woody components.
       ========================================================= */

    /* N flux into new ring (immobile component -> structrual components) */
    s->n_to_alloc_stemimm = s->cstore * f->alstem * p->ncwimm;

    /* N flux into new ring (mobile component -> can be retrans for new
       woody tissue) */
    s->n_to_alloc_stemmob = s->cstore * f->alstem * (p->ncwnew - p->ncwimm);
    s->n_to_alloc_branch = s->cstore * f->albranch * p->ncbnew;
    s->n_to_alloc_croot = s->cstore * f->alcroot * p->nccnew;

    /* Calculate remaining N left to allocate to leaves and roots */
    ntot = MAX(0.0, (s->nstore - s->n_to_alloc_stemimm - s->n_to_alloc_stemmob -
                     s->n_to_alloc_branch));

    /* allocate remaining N to flexible-ratio pools */
    s->n_to_alloc_shoot = (ntot * f->alleaf /
                            (f->alleaf + f->alroot * p->ncrfac));
    s->n_to_alloc_root = ntot - s->n_to_alloc_shoot;

    /*
    leaf_NC = s->n_to_alloc_shoot / s->c_to_alloc_shoot
    if leaf_NC > 0.04:
        s->n_to_alloc_shoot = s->c_to_alloc_shoot * 0.14

    s->n_to_alloc_root = ntot - s->n_to_alloc_shoot


    root_NC = s->n_to_alloc_root / s->c_to_alloc_root
    ncmaxr = 0.04 * p->ncrfac
    if root_NC > ncmaxr:
        extrar = (s->n_to_alloc_root -
                  (s->c_to_alloc_root * ncmaxr))

        s->inorgn += extrar
        s->n_to_alloc_root -= extrar
    */
    
    /* =========================================================
    Phosphorus - Fixed ratios P allocation to woody components.
    ========================================================= */
    
    /* P flux into new ring (immobile component -> structrual components) */
    s->p_to_alloc_stemimm = s->cstore * f->alstem * p->pcwimm;
    
    /* P flux into new ring (mobile component -> can be retrans for new
    woody tissue) */
    s->p_to_alloc_stemmob = s->cstore * f->alstem * (p->pcwnew - p->pcwimm);
    s->p_to_alloc_branch = s->cstore * f->albranch * p->pcbnew;
    s->p_to_alloc_croot = s->cstore * f->alcroot * p->pccnew;
    
    /* Calculate remaining P left to allocate to leaves and roots */
    ptot = MAX(0.0, (s->pstore - s->p_to_alloc_stemimm - s->p_to_alloc_stemmob -
    s->p_to_alloc_branch));
    
    /* allocate remaining P to flexible-ratio pools */
    s->p_to_alloc_shoot = (ptot * f->alleaf /
                             (f->alleaf + f->alroot * p->pcrfac));
    s->p_to_alloc_root = ptot - s->p_to_alloc_shoot;
    
    
    return;
}


double nitrogen_retrans(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
    /* Nitrogen retranslocated from senesced plant matter.
    Constant rate of n translocated from mobile pool

    Parameters:
    -----------
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate

    Returns:
    --------
    N retrans : float
        N retranslocated plant matter

    */
    double leafretransn, rootretransn, crootretransn, branchretransn,
           stemretransn;

    if (c->deciduous_model) {
        leafretransn = p->fretrans * f->lnrate * s->remaining_days[doy];
    } else {
        leafretransn = p->fretrans * fdecay * s->shootn;
    }

    rootretransn = p->rretrans * rdecay * s->rootn;
    crootretransn = p->cretrans * p->crdecay * s->crootn;
    branchretransn = p->bretrans * p->bdecay * s->branchn;
    stemretransn = (p->wretrans * p->wdecay * s->stemnmob + p->retransmob *
                    s->stemnmob);

    /* store for NCEAS output */
    f->leafretransn = leafretransn;

    return (leafretransn + rootretransn + crootretransn + branchretransn +
            stemretransn);
}

double phosphorus_retrans(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
  /* Phosphorus retranslocated from senesced plant matter.
   Constant rate of p translocated from mobile pool
   
   Parameters:
   -----------
   fdecay : float
   foliage decay rate
   rdecay : float
   fine root decay rate
   
   Returns:
   --------
   P retrans : float
   P retranslocated plant matter
   
   */
  double leafretransp, rootretransp, crootretransp, branchretransp,
  stemretransp;
  
  if (c->deciduous_model) {
    leafretransp = p->fretransp * f->lprate * s->remaining_days[doy];
  } else {
    leafretransp = p->fretransp * fdecay * s->shootp;
  }
  
  rootretransp = p->rretrans * rdecay * s->rootp;
  crootretransp = p->cretrans * p->crdecay * s->crootp;
  branchretransp = p->bretrans * p->bdecay * s->branchp;
  stemretransp = (p->wretrans * p->wdecay * s->stempmob + p->retransmob *
    s->stempmob);
  
  /* store for NCEAS output */
  f->leafretransp = leafretransp;
  
  return (leafretransp + rootretransp + crootretransp + branchretransp +
          stemretransp);
}

double calculate_nuptake(control *c, params *p, state *s) {
    /* N uptake depends on the rate at which soil mineral N is made
    available to the plants.

    Returns:
    --------
    nuptake : float
        N uptake

    References:
    -----------
    * Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.
    * Raich et al. 1991, Ecological Applications, 1, 399-429.

    */
    double nuptake, U0, Kr;

    if (c->nuptake_model == 0) {
        /* Constant N uptake */
        nuptake = p->nuptakez;
    } else if (c->nuptake_model == 1) {
        /* evaluate nuptake : proportional to dynamic inorganic N pool */
        nuptake = p->rateuptake * s->inorgn;
    } else if (c->nuptake_model == 2) {
        /* N uptake is a saturating function on root biomass following
           Dewar and McMurtrie, 1996. */

        /* supply rate of available mineral N */
        U0 = p->rateuptake * s->inorgn;
        Kr = p->kr;
        nuptake = MAX(U0 * s->root / (s->root + Kr), 0.0);

        /* Make minimum uptake rate supply rate for deciduous_model cases
           otherwise it is possible when growing from scratch we don't have
           enough root mass to obtain N at the annual time step
           I don't see an obvious better solution?
        if c->deciduous_model:
            nuptake = max(U0 * s->root / (s->root + Kr), U0) */
    } else {
        fprintf(stderr, "Unknown N uptake option\n");
        exit(EXIT_FAILURE);
    }

    return (nuptake);
}


double calculate_puptake(control *c, params *p, state *s) {
  /* P uptake depends on the rate at which soil mineral P is made
  available to the plants.
  
  Returns:
  --------
  puptake : float
  P uptake

  
  */
  double puptake, U0, Kr;
  
  if (c->puptake_model == 0) {
    /* Constant P uptake */
    puptake = p->puptakez;
  } else if (c->puptake_model == 1) {
    /* evaluate puptake : proportional to lab P pool that is available to plant uptake (a function of mineral N) */
    puptake = p->prateuptake * s->inorglabp * p->p_lab_avail;
  } else if (c->puptake_model == 2) {
    /* P uptake is a saturating function on root biomass following
    Dewar and McMurtrie, 1996. */
    
    /* supply rate of available mineral P */
    U0 = p->prateuptake * s->inorglabp * p->p_lab_avail;
    Kr = p->krp;
    puptake = MAX(U0 * s->root / (s->root + Kr), 0.0);
    
    /* Make minimum uptake rate supply rate for deciduous_model cases
    otherwise it is possible when growing from scratch we don't have
    enough root mass to obtain N at the annual time step
    I don't see an obvious better solution?
    if c->deciduous_model:
    nuptake = max(U0 * s->root / (s->root + Kr), U0) */
  } else {
    fprintf(stderr, "Unknown P uptake option\n");
    exit(EXIT_FAILURE);
  }
  
  return (puptake);
}





