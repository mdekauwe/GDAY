#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "gday.h"

typedef struct {
    FILE *ifp;
    FILE *ofp;
    FILE *ofp_sd;
    FILE *ofp_hdr;
    char  cfg_fname[STRING_LENGTH];
    char  met_fname[STRING_LENGTH];
    char  out_fname[STRING_LENGTH];
    char  out_subdaily_fname[STRING_LENGTH];
    char  out_fname_hdr[STRING_LENGTH];
    char  out_param_fname[STRING_LENGTH];
    char  git_hash[STRING_LENGTH];
    int   adjust_rtslow;
    int   alloc_model;
    int   assim_model;
    int   calc_sw_params;
    int   deciduous_model;
    int   disturbance;
    int   fixed_stem_nc;
    int   fixed_stem_pc;
    int   fixed_lai;
    int   fixleafnc;
    int   fixleafpc;
    int   grazing;
    int   gs_model;
    int   model_optroot;
    int   modeljm;
    int   ncycle;
    int   pcycle;
    int   num_years;
    int   nuptake_model;
    int   puptake_model;
    int   output_ascii;
    int   passiveconst;
    int   print_options;
    int   ps_pathway;
    int   respiration_model;
    int   strfloat;
    int   strpfloat;
    int   sw_stress_model;
    int   text_effect_p;
    int   use_eff_nc;
    int   water_stress;
    int   num_days;
    int   total_num_days;
    char  git_code_ver[STRING_LENGTH];
    int   spin_up;
    int   PRINT_GIT;
    int   hurricane;
    int   exudation;
    int   sub_daily;
    int   num_hlf_hrs;
    long  hour_idx;
    long  day_idx;

} control;


typedef struct {
    double activesoil;                  /* active C som pool (t/ha) */
    double activesoiln;                 /* active N som pool (t/ha) */
    double activesoilp;                 /* active P som pool (t/ha) */
    double age;                         /* Current stand age (years) */
    double avg_albranch;                /* Average branch growing season allocation fractions */
    double avg_alcroot;                 /* Average coarse root growing season allocation fractions */
    double avg_alleaf;                  /* Average leaf growing season allocation fractions */
    double avg_alroot;                  /* Average fine root growing season allocation fractions */
    double avg_alstem;                  /* Average stem growing season allocation fractions */
    double branch;                      /* branch c (t/ha) */
    double branchn;                     /* branch n (t/ha) */
    double branchp;                     /* branch p (t/ha) */
    double canht;                       /* canopy height (m) */
    double croot;                       /* coarse root c (t/ha) */
    double crootn;                      /* coarse root n (t/ha) */
    double crootp;                      /* coarse root p (t/ha) */
    double cstore;                      /* C store for deciduous model (t/ha) annual ? */
    double inorgn;                      /* Inorganic soil N pool - dynamic (t/ha) */
    double inorgp;                      /* Inorganic soil P pool - dynamic (t/ha) */
    double inorgminp;                      /* Inorganic soil P pool - mineral P = lab + sorbed (t/ha) */
    double inorglabp;                   /* Inorganic soil P pool - labile P (t/ha) */
    double inorgsorbp;                   /* Inorganic soil P pool - sorbed P (t/ha) */
    double inorgssorbp;                   /* Inorganic soil P pool - strongly sorbed P (t/ha) */
    double inorgoccp;                   /* Inorganic soil P pool - occluded P (t/ha) */
    double inorgparp;                   /* Inorganic soil P pool - parent P (t/ha) */
    double lai;                         /* leaf area index m2 (leaf) m-2 (ground) */
    double fipar;
    double metabsoil;                   /* metabolic soil c (t/ha) */
    double metabsoiln;                  /* metabolic soil n (t/ha) */
    double metabsoilp;                  /* metabolic soil p (t/ha) */
    double metabsurf;                   /* metabolic surface c (t/ha) */
    double metabsurfn;                  /* metabolic surface n (t/ha) */
    double metabsurfp;                  /* metabolic surface p (t/ha) */
    double nstore;                      /* N store for deciduous model (t/ha) */
    double pstore;                      /* P store for deciduous model (t/ha) */
    double passivesoil;                 /* passive C som pool (t/ha) */
    double passivesoiln;                /* passive N som pool (t/ha) */   
    double passivesoilp;                /* passive P som pool (t/ha) */
    double pawater_root;                /* plant available water - root zone (mm) */
    double pawater_topsoil;             /* plant available water - top soil(mm) */
    double prev_sma;
    double root;                        /* root c (t/ha) */
    double root_depth;                  /* rooting depth, Dmax (m) */
    double rootn;                       /* root n (t/ha) */
    double rootp;                       /* root p (t/ha) */
    double sapwood;
    double shoot;                       /* shoot c (t/ha) */
    double shootn;                      /* shoot n (t/ha) */
    double shootp;                      /* shoot p (t/ha) */
    double sla;                         /* specific leaf area */
    double slowsoil;                    /* slow C som pool (t/ha) */
    double slowsoiln;                   /* slow N som pool (t/ha) */
    double slowsoilp;                   /* slow P som pool (t/ha) */
    double stem;
    double stemn;                       /* Stem N (t/ha) = stemnimm + stemnmob */
    double stemnimm;
    double stemnmob;
    double stemp;                       /* Stem P (t/ha) = stempimm + stempmob */
    double stempimm;
    double stempmob;
    double structsoil;                  /* soil structural c (t/ha) */
    double structsoiln;                 /* soil structural n (t/ha) */
    double structsoilp;                 /* soil structural p (t/ha) */
    double structsurf;                  /* surface structural c (t/ha) */
    double structsurfn;                 /* surface structural n (t/ha) */
    double structsurfp;                 /* surface structural p (t/ha) */
    double shootnc;                     /* shoot nc ratio */
    double rootnc;                      /* root pn ratio */
    double shootpc;                     /* shoot pc ratio */
    double rootpc;                      /* root pc ratio */
    double remaining_days[366];
    double wtfac_root;
    double wtfac_topsoil;
    double delta_sw_store;
    double b_root;
    double theta_sat_root;
    double theta_sat_topsoil;
    double psi_sat_root;
    double psi_sat_topsoil;
    double b_topsoil;
    double leaf_out_days[366];
    double growing_days[366];
    double c_to_alloc_shoot;
    double n_to_alloc_shoot;
    double p_to_alloc_shoot;
    double c_to_alloc_root;
    double n_to_alloc_root;
    double p_to_alloc_root;
    double c_to_alloc_croot;
    double n_to_alloc_croot;
    double p_to_alloc_croot;
    double c_to_alloc_branch;
    double n_to_alloc_branch;
    double p_to_alloc_branch;
    double c_to_alloc_stem;
    double n_to_alloc_stemmob;
    double n_to_alloc_stemimm;
    double p_to_alloc_stemmob;
    double p_to_alloc_stemimm;
    double anpp;                    /* aboveground NPP */
    double litterc;                 /* litter carbon */
    double littern;                 /* litter nitrogen */
    double litterp;                 /* litter phosphorus */
    double littercbg;               /* litter C belowground */
    double littercag;               /* litter C aboveground */
    double litternag;               /* litter N aboveground */
    double litternbg;               /* litter N belowground */
    double litterpag;               /* litter P aboveground */
    double litterpbg;               /* litter P belowground */
    double plantc;                  /* plant C */
    double plantn;                  /* plant N */
    double plantp;                  /* plant P */
    double totalc;                  /* total C */
    double totaln;                  /* total N */
    double totalp;                  /* total P */
    double soilc;                   /* Soil C */
    double soiln;                   /* soil N */
    double soilp;                   /* soil P */
    double canopy_store;
    double psi_s_topsoil;
    double psi_s_root;
} state;

typedef struct {
    double a0rhizo;                         /* minimum allocation to rhizodeposition [0.0-0.1] */
    double a1rhizo;                         /* slope of allocation to rhizodeposition [0.2-1] */
    double actncmax;                        /* Active pool (=1/3) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double actncmin;                        /* Active pool (=1/15) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double actpcmax;                        /* Active pool (=1/30) P:C ratio of new SOM - maximum [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993) */
    double actpcmin;                        /* Active pool (=1/80) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC]. Based on forest version of CENTURY (Parton et al. 1993) */
    double adapt;
    double ageold;                          /* Plant age when max leaf N C ratio is lowest */
    double ageyoung;                        /* Plant age when max leaf N C ratio is highest */
    double albedo;
    double alpha_c4;                        /* quantium efficiency for C4 plants has no Ci and temp dependancy, so if a fixed constant. */
    double alpha_j;                         /* initial slope of rate of electron transport, used in calculation of quantum yield. Value calculated by Belinda */
    double b_root;
    double b_topsoil;
    double bdecay;                          /* branch and large root turnover rate (1/yr) */
    double biochemical_p_constant;          /* Michaelis-Menton constant for biochemical P mineralisation [g N (g P)-1]; Wang et al., 2007, GB1018*/
    double branch0;                         /* constant in branch-stem allometry (trees) */
    double branch1;                         /* exponent in branch-stem allometry */
    double bretrans;                        /* branch n retranslocation fraction */
    double c_alloc_bmax;                    /* allocation to branches at branch n_crit and p_crit. If using allometric model this is the max alloc to branches */
    double c_alloc_bmin;                    /* allocation to branches at zero branch n/c and p/c. If using allometric model this is the min alloc to branches */
    double c_alloc_cmax;                    /* allocation to coarse roots at n_crit and p_crit. If using allometric model this is the max alloc to coarse roots */
    double c_alloc_fmax;                    /* allocation to leaves at leaf n_crit and p_crit. If using allometric model this is the max alloc to leaves */
    double c_alloc_fmin;                    /* allocation to leaves at zero leaf n/c and p/c. If using allometric model this is the min alloc to leaves */
    double c_alloc_rmax;                    /* allocation to roots at root n_crit and p_crit. If using allometric model this is the max alloc to fine roots */
    double c_alloc_rmin;                    /* allocation to roots at zero root n/c and p/c. If using allometric model this is the min alloc to fine roots */
    double cfracts;                         /* carbon fraction of dry biomass */
    double crdecay;                         /* coarse roots turnover rate (1/yr) */
    double cretrans;                        /* coarse root n retranslocation fraction */
    double croot0;                          /* constant in coarse_root-stem allometry (trees) */
    double croot1;                          /* exponent in coarse_root-stem allometry */
    double crit_n_cost_of_p;                /* Critical value of N cost of root P uptake above which phosphatase production starts [g N (g P)-1]; Wang et al., 2007, GB1018*/
    double ctheta_root;                     /* Fitted parameter based on Landsberg and Waring */
    double ctheta_topsoil;                  /* Fitted parameter based on Landsberg and Waring */
    double cue;                             /* carbon use efficiency, or the ratio of NPP to GPP */
    double d0;
    double d0x;                             /* Length scale for exponential decline of Umax(z) */
    double decayrate[7];
    double d1;
    double delsj;                           /* Deactivation energy for electron transport (J mol-1 k-1) */
    double density;                         /* sapwood density kg DM m-3 (trees) */
    double direct_frac;                     /* direct beam fraction of incident radiation - this is only used with the BEWDY model */
    double displace_ratio;                  /* Value for coniferous forest (0.78) from Jarvis et al 1976, taken from Jones 1992 pg 67. More standard assumption is 2/3 */
    int    disturbance_doy;
    double dz0v_dh;                         /* Rate of change of vegetation roughness length for momentum with height. Value from Jarvis? for conifer 0.075 */
    double eac;                             /* Activation energy for carboxylation [J mol-1] */
    double eag;                             /* Activation energy at CO2 compensation point [J mol-1] */
    double eaj;                             /* Activation energy for electron transport (J mol-1) */
    double eao;                             /* Activation energy for oxygenation [J mol-1] */
    double eav;                             /* Activation energy for Rubisco (J mol-1) */
    double edj;                             /* Deactivation energy for electron transport (J mol-1) */
    double faecescn;                        /* Faeces C:N ratio */
    double faecesn;                         /* Faeces N content */
    double faecescp;                        /* Faeces C:P ratio */
    double faecesp;                         /* Faeces P content */
    double fdecay;                          /* foliage turnover rate (1/yr) */
    double fdecaydry;                       /* Foliage turnover rate - dry soil (1/yr) */
    double fhw;                             /* n:c ratio of stemwood - immobile pool and new ring */
    double fhwp;                            /* p:c ratio of stemwood - immobile pool and new ring */
    double finesoil;                        /* clay+silt fraction */
    double fmleaf;
    double fmroot;
    double fmfaeces;
    double fix_lai;                         /* value to fix LAI to, control fixed_lai flag must be set */
    double fracfaeces;                      /* Fractn of grazd C that ends up in faeces (0..1) */
    double fracteaten;                      /* Fractn of leaf production eaten by grazers */
    double fractosoil;                      /* Fractn of grazed N recycled to soil:faeces+urine */
    double fractosoilp;                     /* Fractn of grazed P recycled to soil:faeces+urine */
    double fractup_soil;                    /* fraction of uptake from top soil layer */
    double fretrans;                        /* foliage n retranslocation fraction - 46-57% in young E. globulus trees - see Corbeels et al 2005 ecological modelling 187, pg 463. Roughly 50% from review Aerts '96 */
    double fretransp;                       /* foliage p retranslocation fraction - 39.5-69 in Southern US FACE site - Finzi et al. 2001 Ecology  */
    double g1;                              /* stomatal conductance parameter: Slope of reln btw gs and assimilation (fitted by species/pft). */
    double gamstar25;                       /* Base rate of CO2 compensation point at 25 deg C [umol mol-1] */
    double growth_efficiency;               /* growth efficiency (yg) - used only in Bewdy */
    double height0;                         /* Height when leaf:sap area ratio = leafsap0 (trees) */
    double height1;                         /* Height when leaf:sap area ratio = leafsap1 (trees) */
    double heighto;                         /* constant in avg tree height (m) - stem (t C/ha) reln */
    double htpower;                         /* Exponent in avg tree height (m) - stem (t C/ha) reln */
    double intercep_frac;                   /* Maximum intercepted fraction, values in Oishi et al 2008, AFM, 148, 1719-1732 ~13.9% +/- 4.1, so going to assume 15% following Landsberg and Sands 2011, pg. 193. */
    double jmax;                            /* maximum rate of electron transport (umol m-2 s-1) */
    double jmaxna;                          /* slope of the reln btween jmax and leaf N content, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010. */
    double jmaxnb;                          /* intercept of jmax vs n, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010. */
    double jv_intercept;                    /* Jmax to Vcmax intercept */
    double jv_slope;                        /* Jmax to Vcmax slope */
    double kc25;                            /* Base rate for carboxylation by Rubisco at 25degC [mmol mol-1] */
    double kdec1;                           /* surface structural decay rate (1/yr) */
    double kdec2;                           /* surface metabolic decay rate (1/yr) */
    double kdec3;                           /* soil structural decay rate (1/yr) */
    double kdec4;                           /* soil metabolic decay rate(1/yr) */
    double kdec5;                           /* active pool decay rate (1/yr) */
    double kdec6;                           /* slow pool decay rate (1/yr) */
    double kdec7;                           /* passive pool decay rate (1/yr) */
    double kext;                            /* extinction coefficient for function calculate_top_of_canopy_n/p in photosynthesis.c */
    double kn;                              /* extinction coefficient of nitrogen in the canopy, assumed to be 0.3 by defaul which comes half from Belinda's head and is supported by fig 10 in Lloyd et al. Biogeosciences, 7, 1833–1859, 2010 */
    double kp;                              /* extinction coefficient of phosphorus in the canopy for function calculate_top_of_canopy_leafp in canopy.c */
    double knl;
    double ko25;                            /* Base rate for oxygenation by Rubisco at 25degC [umol mol-1]. Note value in Bernacchie 2001 is in mmol!! */
    double kq10;                            /* exponential coefficient for Rm vs T */
    double kr;                              /* N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha) see Silvia's PhD, Dewar and McM, 96. */
    double krp;                             /* P uptake coefficent */
    double ks;                              /* an empirical constant [t P ha-1] - sorption capacity increased with the age of the substrate */
    double lad;                             /* Leaf angle distribution: 0 = spherical leaf angle distribution; 1 = horizontal leaves; -1 = vertical leaves */
    double lai_closed;                      /* LAI of closed canopy (max cover fraction is reached (m2 (leaf) m-2 (ground) ~ 2.5) */
    double latitude;                        /* latitude (degrees, negative for south) */
    double longitude;                       /* longitude (degrees, negative for west) */
    double leafsap0;                        /* leaf area  to sapwood cross sectional area ratio when Height = Height0 (mm^2/mm^2) */
    double leafsap1;                        /* leaf to sap area ratio when Height = Height1 (mm^2/mm^2) */
    double ligfaeces;                       /* Faeces lignin as fraction of biomass */
    double ligroot;                         /* lignin-to-biomass ratio in root litter; Values from White et al. = 0.22  - Value in Smith et al. 2013 = 0.16, note subtly difference in eqn C9. */
    double ligshoot;                        /* lignin-to-biomass ratio in leaf litter; Values from White et al. DBF = 0.18; ENF = 0.24l; GRASS = 0.09; Shrub = 0.15 - Value in smith et al. 2013 = 0.2, note subtly difference in eqn C9. */
    double liteffnc;
    double max_intercep_lai;                /* canopy LAI at which interception is maximised. */
    double max_p_biochemical;               /* max rate of biochemical P mineralisation [g P m-2 y-1]; Wang et al., 2007, GB1018*/
    double measurement_temp;                /* temperature Vcmax/Jmax are measured at, typical 25.0 (celsius)  */
    double ncbnew;                          /* N alloc param: new branch N C at critical leaf N C */
    double ncbnewz;                         /* N alloc param: new branch N C at zero leaf N C */
    double nccnew;                          /* N alloc param: new coarse root N C at critical leaf N C */
    double nccnewz;                         /* N alloc param: new coarse root N C at zero leaf N C */
    double ncmaxfold;                       /* max N:C ratio of foliage in old stand, if the same as young=no effect */
    double ncmaxfyoung;                     /* max N:C ratio of foliage in young stand, if the same as old=no effect */
    double ncmaxr;                          /* max N:C ratio of roots */
    double ncrfac;                          /* N:C of fine root prodn / N:C of leaf prodn */
    double ncwimm;                          /* N alloc param: Immobile stem N C at critical leaf N C */
    double ncwimmz;                         /* N alloc param: Immobile stem N C at zero leaf N C */
    double ncwnew;                          /* N alloc param: New stem ring N:C at critical leaf N:C (mob) */
    double ncwnewz;                         /* N alloc param: New stem ring N:C at zero leaf N:C (mobile) */
    double nf_crit;                         /* leaf N:C below which N availability limits productivity  */
    double nf_min;                          /* leaf N:C minimum N concentration which allows productivity */
    double nmax;
    double nmin;                            /* (bewdy) minimum leaf n for +ve p/s (g/m2) */
    double nmin0;                           /* mineral N pool corresponding to Actnc0,etc (g/m2) */
    double nmincrit;                        /* Critical mineral N pool at max soil N:C (g/m2) (Parton et al 1993, McMurtrie et al 2001). */
    double ntheta_root;                     /* Fitted parameter based on Landsberg and Waring */
    double ntheta_topsoil;                  /* Fitted parameter based on Landsberg and Waring */
    double nuptakez;                        /* constant N uptake per year (1/yr) */
    double oi;                              /* intercellular concentration of O2 [umol mol-1] */
    double p_atm_deposition;                /* atmospheric P deposition rate [t/ha/yr] Newman 1995, Journal of Ecology Table 5 (Tennessee) */
    double p_rate_par_weather;              /* parent P material weathering rate [yr-1] */
    double passivesoilnz;
    double passivesoilpz;
    double passivesoilz;
    double passncmax;                       /* Passive pool (=1/7) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double passncmin;                       /* Passive pool (=1/10) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double passpcmax;                       /* Passive pool (=1/20) P:C ratio of new SOM - maximum [units: gP/gC] */
    double passpcmin;                       /* Passive pool (=1/200) P:C of new SOM - when Pmin=Pmin0 [units: gP/gC] */
    double pcbnew;                          /* P alloc param: new branch P C at critical leaf P C */
    double pcbnewz;                         /* P alloc param: new branch P C at zero leaf P C */
    double pccnew;                          /* P alloc param: new coarse root P C at critical leaf P C */
    double pccnewz;                         /* P alloc param: new coarse root P C at zero leaf P C */
    double pcmaxfold;                       /* max P:C ratio of foliage in old stand, if the same as young=no effect */
    double pcmaxfyoung;                     /* max P:C ratio of foliage in young stand, if the same as old=no effect */
    double pcmaxr;                          /* max P:C ratio of roots */
    double pcrfac;                          /* P:C of fine root prodp / P:C c of leaf prodp */
    double pcwimm;                          /* P alloc param: Immobile stem P C at critical leaf P C */
    double pcwimmz;                         /* P alloc param: Immobile stem P C at zero leaf P C */
    double pcwnew;                          /* P alloc param: New stem ring P:C at critical leaf P:C (mob) */
    double pcwnewz;                         /* P alloc param: New stem ring P:C at zero leaf P:C (mobile) */
    double pf_crit;                         /* leaf P:C below which P availability limits productivity  */
    double pf_min;                          /* leaf P:C minimum P concentration which allows productivity */
    double phmax;                           /* max pH for determining effect on solubility of secondary P */
    double phmin;                           /* min pH for determining effect on solubility of secondary P */
    double phtextmin;                       /* the solubility of secondary P corresponding to min pH (/yr) */
    double phtextmax;                       /* the solubility of secondary P corresponding to max pH (/yr) */
    double phtextslope;                     /* slope controlling effect of sand on secondary P flow to mineral P */
    double p_lab_avail;                     /* Fraction of labile P available for plant uptake */
    double pmax;
    double pmin;                            /* (bewdy) minimum leaf p for +ve p/s (g/m2) */
    double pmin0;                           /* mineral P pool corresponding to Actpc0,etc (g/m2) */
    double pmincrit;                        /* Critical mineral P pool at max soil P:C (g/m2) */
    double prateloss;                       /* Rate of P loss from mineral P pool (/yr), Ref Wang et al., 2007, GB1018 */
    double prateuptake;                     /* Rate of P uptake from mineral P pool (/yr), guess value */
    double prescribed_leaf_NC;              /* If the N-Cycle is switched off this needs to be set, e.g. 0.03 */
    double prescribed_leaf_PC;              /* If the P-Cycle is switched off this needs to be set, e.g. 0.00249 */
    double previous_ncd;                    /* In the first year we don't have last years data, so I have precalculated the average of all the november-jan chilling values  */
    double psecmnp;                         /* controls the flow from secondary to mineral P, used when text_effect_p set to 0 */
    double psi_sat_root;                    /* MPa */
    double psi_sat_topsoil;                 /* MPa */
    double psie_topsoil;                    /* Soil water potential at saturation (m) */
    double psie_root;                       /* Soil water potential at saturation (m) */
    double puptakez;                        /* constant P uptake per year (1/yr) */
    double qs;                              /* exponent in water stress modifier, =1.0 JULES type representation, the smaller the values the more curved the depletion.  */
    double r0;                              /* root C at half-maximum N uptake (kg C/m3) */
    double rate_ssorb_occ;                  /* Rate constant of the transfer of P from strongly sorbed pool to occluded pool, m-1 Yang et al. 2014, Biogeosciences */
    double rate_sorb_ssorb;                 /* Rate constant of the transfer of P from sorbed pool to strongly sorbed pool, m-1 Yang et al. 2014, Biogeosciences */
    double rateloss;                        /* Rate of N loss from mineral N pool (/yr) */
    double rateuptake;                      /* Rate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf Seems to correspond to very low NPP values */
    double rdecay;                          /* root turnover rate (1/yr) */
    double rdecaydry;                       /* root turnover rate - dry soil (1/yr) */
    double retransmob;                      /* Fraction stem mobile N retranscd (/yr) */
    double rfmult;
    double rooting_depth;                   /* Rooting depth (mm) */
    char   rootsoil_type[STRING_LENGTH];
    char   soil_order[STRING_LENGTH];       /* soil order */
    double rretrans;                        /* root n retranslocation fraction */
    //double sand_frac;                       /* fraction of sand in soil (top + root averaged) */ 
    double sapturnover;                     /* Sapwood turnover rate: conversion of sapwood to heartwood (1/yr) */
    double sla;                             /* specific leaf area (m2 one-sided/kg DW) */
    double slamax;                          /* (if equal slazero=no effect) specific leaf area new fol at max leaf N/C (m2 one-sided/kg DW) */
    double slazero;                         /* (if equal slamax=no effect) specific leaf area new fol at zero leaf N/C (m2 one-sided/kg DW) */
    double slowncmax;                       /* Slow pool (=1/15) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double slowncmin;                       /* Slow pool (=1/40) N:C of new SOM - when Nmin=Nmin0" [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double slowpcmax;                       /* Slow pool (=1/90) P:C ratio of new SOM - maximum [units: gP/gC]. */
    double slowpcmin;                       /* Slow pool (=1/200) P:C of new SOM - when Pmin=Pmin0" [units: gP/gC]. */
    double smax;                            /* Maximum amount of sorbed P in the soil [tt P ha-1] */
    double store_transfer_len;
    double structcn;                        /* C:N ratio of structural bit of litter input */
    double structrat;                       /* structural input n:c as fraction of metab */
    double structcp;                        /* C:P ratio of structural bit of litter input, Ref Attiwill 1980, Aus. J. Bot. 28, 199-222 Table 9 sum of branch, stem, sap and heartwood; */
    double structratp;                      /* structural input p:c as fraction of metab */
    double soilph;                          /* soil pH value */
    double sorpmx;                          /* maximum P sorption potential for a soil */
    double sorpaf;                          /* slope term which controls the fraction of mineral P that is labile */
    double targ_sens;                       /* sensitivity of allocation (leaf/branch) to track the target, higher values = less responsive. */
    double theta;                           /* curvature of photosynthetic light response curve */
    double theta_sp_root;
    double theta_sp_topsoil;
    double theta_fc_root;
    double theta_fc_topsoil;
    double theta_wp_root;
    double theta_wp_topsoil;
    double topsoil_depth;                   /* Topsoil depth (mm) */
    char   topsoil_type[STRING_LENGTH];
    double vcmax;                           /* maximum rate of carboxylation (umol m-2 s-1)  */
    double vcmaxna;                         /* slope of the reln btween vcmax and leaf N content, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010. */
    double vcmaxnb;                         /* intercept of vcmax vs n, units = (umol [gN]-1 s-1) # And for Vcmax-N slopes (vcmaxna) see Table 8.2 in CLM4_tech_note, Oleson et al. 2010. */
    double watdecaydry;                     /* water fractn for dry litterfall rates */
    double watdecaywet;                     /* water fractn for wet litterfall rates */
    double wcapac_root;                     /* Max plant avail soil water -root zone, i.e. total (mm) (smc_sat-smc_wilt) * root_depth (750mm) = [mm (water) / m (soil depth)] */
    double wcapac_topsoil;                  /* Max plant avail soil water -top soil (mm) */
    double wdecay;                          /* wood turnover rate (1/yr) */
    double wetloss;                         /* Daily rainfall lost per lai (mm/day) */
    double wretrans;                        /* mobile wood N retranslocation fraction */
    double z0h_z0m;                         /* Assume z0m = z0h, probably a big assumption [as z0h often < z0m.], see comment in code!! But 0.1 might be a better assumption */
    int    growing_seas_len;
    double prime_y;
    double prime_z;
    int    return_interval;                 /* years */
    int    burn_specific_yr;
    int    hurricane_doy;
    int    hurricane_yr;
    double root_exu_CUE;
    double leaf_width;
    double leaf_abs;
} params;

typedef struct {

    double *year;
    double *rain;
    double *par;
    double *tair;
    double *tsoil;
    double *co2;
    double *ndep;
    double *pdep;
    double *nfix;       /* N inputs from biological fixation (t/ha/timestep (d/30min)) */
    double *pfix;       /* P inputs from biological fixation (t/ha/timestep (d/30min)) */
    double *wind;
    double *press;

    /* Day timestep */
    double *prjday; /* should really be renamed to doy for consistancy */
    double *tam;
    double *tpm;
    double *tmin;
    double *tmax;
    double *tday;
    double *vpd_am;
    double *vpd_pm;
    double *wind_am;
    double *wind_pm;
    double *par_am;
    double *par_pm;

    /* sub-daily timestep */
    double *vpd;
    double *doy;
    double *diffuse_frac;


} met_arrays;


typedef struct {

    /* sub-daily */
    double rain;
    double wind;
    double press;
    double vpd;
    double tair;
    double sw_rad;
    double par;
    double Ca;
    double ndep;
    double pdep;
    double nfix;       /* N inputs from biological fixation (t/ha/timestep (d/30min)) */
    double Pfix;       /* P inputs from biological fixation (t/ha/timestep (d/30min)) */
    double tsoil;

    /* daily */
    double tair_am;
    double tair_pm;
    double sw_rad_am;
    double sw_rad_pm;
    double vpd_am;
    double vpd_pm;
    double wind_am;
    double wind_pm;
    double Tk_am;
    double Tk_pm;

} met;



typedef struct {
    double gpp_gCm2;
    double npp_gCm2;
    double gpp_am;
    double gpp_pm;
    double gpp;
    double npp;
    double nep;
    double auto_resp;
    double hetero_resp;
    double retrans;         /* plnat n retranslocation */
    double retransp;        /* plant p retranslocation */
    double apar;

    /* n */
    double nuptake;         /* n plant uptake rate */
    double nloss;           /* n loss by leaching and volatilisation */
    double npassive;        /* n passive -> active */
    double ngross;          /* N gross mineralisation */
    double nimmob;          /* N immobilisation in SOM */
    double nlittrelease;    /* N rel litter = struct + metab */
    double activelossf;     /* frac of active C -> CO2 */
    double nmineralisation;
    
    /* p */
    double puptake;         /* P plant uptake rate */
    double ploss;           /* P loss by leaching */
    double ppassive;        /* P passive -> active */
    double pgross;          /* P gross mineralisation */
    double pimmob;          /* P immobilisation in SOM */
    double plittrelease;    /* P rel litter = struct + metab */
    double pmineralisation; /* P net mineralised */
    

    /* water fluxes */
    double wue;
    double et;
    double soil_evap;
    double transpiration;
    double interception;
    double throughfall;
    double canopy_evap;
    double runoff;
    double gs_mol_m2_sec;
    double ga_mol_m2_sec;
    double omega;

    /* daily C production */
    double cpleaf;
    double cproot;
    double cpcroot;
    double cpbranch;
    double cpstem;

    /* daily N production */
    double npleaf;
    double nproot;
    double npcroot;
    double npbranch;
    double npstemimm;
    double npstemmob;
    double nrootexudate;
    
    /* daily P production */
    double ppleaf;
    double pproot;
    double ppcroot;
    double ppbranch;
    double ppstemimm;
    double ppstemmob;
    double prootexudate;

    /* dying stuff */
    double deadleaves;      /* Leaf litter C production (t/ha/yr) */
    double deadroots;       /* Root litter C production (t/ha/yr) */
    double deadcroots;      /* Coarse root litter C production (t/ha/yr) */
    double deadbranch;      /* Branch litter C production (t/ha/yr) */
    double deadstems;       /* Stem litter C production (t/ha/yr) */
    double deadleafn;       /* Leaf litter N production (t/ha/yr) */
    double deadrootn;       /* Root litter N production (t/ha/yr) */
    double deadcrootn;      /* Coarse root litter N production (t/ha/yr) */
    double deadbranchn;     /* Branch litter N production (t/ha/yr) */
    double deadstemn;       /* Stem litter N production (t/ha/yr) */
    double deadleafp;       /* Leaf litter P production (t/ha/yr) */
    double deadrootp;       /* Root litter P production (t/ha/yr) */
    double deadcrootp;      /* Coarse root litter P production (t/ha/yr) */
    double deadbranchp;     /* Branch litter P production (t/ha/yr) */
    double deadstemp;       /* Stem litter P production (t/ha/yr) */
    double deadsapwood;

    /* grazing stuff */
    double ceaten;          /* C consumed by grazers (t C/ha/y) */
    double neaten;          /* N consumed by grazers (t N/ha/y) */
    double peaten;          /* P consumed by grazers (t P/ha/y) */
    
    double faecesc;         /* Flux determined by faeces C:N */
    double nurine;          /* Rate of N input to soil in urine (t/ha/y) */
    double purine;          /* Rate of P input to soil in urine (t/ha/y) */
    
    double leafretransn;    /* N retranslocation leaf */
    double leafretransp;    /* P version of leafretransn */
    

    /* C, N  & P Surface litter */
    double surf_struct_litter;
    double surf_metab_litter;
    double n_surf_struct_litter;
    double n_surf_metab_litter;
    double p_surf_struct_litter;
    double p_surf_metab_litter;

    /* C, N & P Root Litter */
    double soil_struct_litter;
    double soil_metab_litter;
    double n_soil_struct_litter;
    double n_soil_metab_litter;
    double p_soil_struct_litter;
    double p_soil_metab_litter;
    

    /* C, N & P litter fluxes to slow pool */
    double surf_struct_to_slow;
    double soil_struct_to_slow;
    double n_surf_struct_to_slow;
    double n_soil_struct_to_slow;
    double p_surf_struct_to_slow;
    double p_soil_struct_to_slow;
    
    /* C, N & P litter fluxes to active pool */
    double surf_struct_to_active;
    double soil_struct_to_active;
    double n_surf_struct_to_active;
    double n_soil_struct_to_active;
    double p_surf_struct_to_active;
    double p_soil_struct_to_active;

    /* Metabolic fluxes to active pool */
    double surf_metab_to_active;
    double soil_metab_to_active;
    double n_surf_metab_to_active;
    double n_soil_metab_to_active;
    double p_surf_metab_to_active;
    double p_soil_metab_to_active;

    /* fluxes out of active pool */
    double active_to_slow;
    double active_to_passive;
    double n_active_to_slow;
    double n_active_to_passive;
    double p_active_to_slow;
    double p_active_to_passive;

    /* Fluxes from slow pool */
    double slow_to_active;
    double slow_to_passive;
    double n_slow_to_active;
    double n_slow_to_passive;
    double p_slow_to_active;
    double p_slow_to_passive;
    double p_slow_biochemical;

    /* fluxes from passive pool */
    double passive_to_active;
    double n_passive_to_active;
    double p_passive_to_active;

    /* C, N & P source fluxes from the active, slow and passive pools */
    double c_into_active;
    double c_into_slow;
    double c_into_passive;
    
    /* inorganic P flux exchanges */
    double p_lab_in;
    double p_lab_out;
    double p_sorb_in;
    double p_sorb_out;
    double p_min_to_ssorb;
    double p_ssorb_to_min;
    double p_ssorb_to_occ;
    double p_par_to_min;
    double p_atm_dep;
    

    /* CO2 flows to the air */
    double co2_to_air[7];

    /* C allocated fracs  */
    double alleaf;             /* allocation to leaf */
    double alroot;             /* allocation to fine root */
    double alcroot;            /* allocation to coarse root */
    double albranch;           /* allocation to branch */
    double alstem;             /* allocation to stems */

    /* Misc stuff */
    double cica_avg; /* used in water balance, only when running mate model */

    double rabove;
    double tfac_soil_decomp;
    double co2_rel_from_surf_struct_litter;
    double co2_rel_from_soil_struct_litter;
    double co2_rel_from_surf_metab_litter;
    double co2_rel_from_soil_metab_litter;
    double co2_rel_from_active_pool;
    double co2_rel_from_slow_pool;
    double co2_rel_from_passive_pool;

    double lrate;
    double wrate;
    double brate;
    double crate;

    /* N allocation rates across growing season */
    double lnrate;
    double bnrate;
    double wnimrate;
    double wnmobrate;
    double cnrate;

    /* P allocation rates across growing season */
    double lprate;              /* leaf allocation rate across growing season */
    double bprate;              /* branch */
    double wpimrate;            /* immobilised */
    double wpmobrate;           /* mobilised */
    double cprate;              /* coarse root */
    
    /* priming/exudation */
    double root_exc;
    double root_exn;
    double root_exp;
    double co2_released_exud;
    double factive;
    double rtslow;
    double rexc_cue;

    double ninflow;  /* N inflow flux from N-fix and Ndep */

} fluxes;

typedef struct {
    /* 2 member arrays are for the sunlit (0) and shaded (1) components */
    int    ileaf;           /* sunlit (0) or shaded (1) leaf index */
    double an_leaf[2];      /* leaf net photosynthesis (umol m-2 s-1) */
    double rd_leaf[2];      /* leaf respiration in the light (umol m-2 s-1) */
    double gsc_leaf[2];     /* leaf stomatal conductance to CO2 (mol m-2 s-1) */
    double apar_leaf[2];    /* leaf abs photosyn. active rad. (umol m-2 s-1) */
    double trans_leaf[2];   /* leaf transpiration (mol m-2 s-1) */
    double rnet_leaf[2];    /* leaf net radiation (W m-2) */
    double lai_leaf[2];     /* sunlit and shaded leaf area (m2 m-2) */
    double omega_leaf[2];   /* leaf decoupling coefficient (-) */
    double tleaf[2];        /* leaf temperature (deg C) */
    double an_canopy;       /* canopy net photosynthesis (umol m-2 s-1) */
    double rd_canopy;       /* canopy respiration in the light (umol m-2 s-1) */
    double gsc_canopy;      /* canopy stomatal conductance to CO2 (mol m-2 s-1) */
    double apar_canopy;     /* canopy abs photosyn. active rad. (umol m-2 s-1) */
    double omega_canopy;    /* canopy decoupling coefficient (-) */
    double trans_canopy;    /* canopy transpiration (mm 30min-1) */
    double rnet_canopy;     /* canopy net radiation (W m-2) */
    double N0;              /* top of canopy nitrogen (g N m-2)) */
    double P0;              /* top of canopy phosphorus (g P m-2)) */
    double elevation;       /* sun elevation angle in degrees */
    double cos_zenith;      /* cos(zenith angle of sun) in radians */
    double diffuse_frac;    /* Fraction of incident rad which is diffuse (-) */
    double direct_frac;     /* Fraction of incident rad which is beam (-) */
    double tleaf_new;       /* new leaf temperature (deg C) */
    double dleaf;           /* leaf VPD (Pa) */
    double Cs;              /* CO2 conc at the leaf surface (umol mol-1) */
    double kb;              /* beam radiation ext coeff of canopy */
    double cscalar[2];      /* scale from single leaf to canopy */
} canopy_wk;


#endif
