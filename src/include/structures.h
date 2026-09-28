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
    int   fixed_lai;
    int   prescribed_lai;  /* LAI read from the met file (extra last column) */
    int   met_start_year;  /* netCDF forcing: first year to use (0 = all) */
    int   met_end_year;    /* netCDF forcing: last year to use (0 = all) */
    int   lai_pft_index;   /* netCDF LAI file: index on the pft dimension */
    char  lai_fname[STRING_LENGTH];  /* netCDF LAI file (optional) */
    char  lai_var[STRING_LENGTH];    /* its LAI variable name */
    int   fixleafnc;
    int   grazing;
    int   gs_model;
    int   model_optroot;
    int   modeljm;
    int   ncycle;
    int   pcycle;   /* phosphorus cycle on or off */
    int   strpfloat;   /* structural litter P:C floats (else structcp) */
    int   puptake_model;   /* P uptake: 0 constant, 1 labile P x rate, 2 as 1 x root/(root+krp) */
    int   text_effect_p;   /* strongly sorbed -> mineral P: 1 CENTURY pH/texture, 0 constant psecmnp */
    int   num_years;
    int   nuptake_model;
    int   output_ascii;
    int   passiveconst;
    int   print_options;
    int   ps_pathway;
    int   respiration_model;
    int   strfloat;
    int   sw_stress_model;
    int   use_eff_nc;
    int   water_stress;
    int   nonstomatal_limitation; /* also scale Vcmax/Jmax by beta */
    int   water_balance;
    int   water_store;       /* retired with Emax, must be false */
    int   gs_opt_search;     /* GS_OPT_FLAT or GS_OPT_GOLDEN */
    int   gs_opt_e;          /* E the optimiser costs: GS_OPT_E_TOTAL (gs in
                                series with the leaf boundary layer, the E
                                the energy balance delivers) or
                                GS_OPT_E_STOMATAL (gs alone, as JULES) */
    int   plant_segments;    /* 1 = single plant conductance, 3 = root, stem
                                and leaf segments in series */
    int   root_depth_model;  /* ROOT_DEPTH_DYNAMIC or ROOT_DEPTH_FIXED */
    long  n_supply_bind;     /* gs_opt leaf steps limited by the soil
                                supply (E <= E_supply), reported at the end */
    int   canopy_ga_model;   /* CANOPY_GA_CABLE or CANOPY_GA_SIMPLE, the
                                canopy air space's conductance to the
                                reference height (sub-daily) */
    int   net_lw_model;      /* NET_LW_LWDOWN or NET_LW_MONTEITH (net
                                radiation of the daily model & sub-daily
                                soil) */
    int   canopy_air_space;  /* sub-daily: leaves see canopy air, coupled to
                                the reference height through the canopy
                                aerodynamic conductance (as CABLE) */
    int   root_radial_resistance; /* include the root radial resistance in
                                     the layer uptake weights (SPA) */
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
    int   pdebug;
    int   spinup_method;
    int   soil_drainage;
    int   soil_hydraulics;   /* SAXTON, VAN_GENUCHTEN or BROOKS_COREY */
    int   soil_evap_model;   /* SOIL_EVAP_GDAY, _JULES or _OR (Decker/CABLE) */
    int   canopy_evap_model; /* CANOPY_EVAP_GDAY or _JULES (sub-daily) */
    int   bound_soil_psi;    /* JULES l_bound_soil_wp: psi_close <= psi <= psi_open */
    int   dry_soil_correction; /* JULES l_ds_correction (Webb 2000) below psi_close */
} control;


typedef struct {
    double *day_length;

    double activesoil;                  /* active C som pool (t/ha) */
    double activesoiln;                 /* active N som pool (t/ha) */
    double age;                         /* Current stand age (years) */
    double avg_albranch;                /* Average branch growing season allocation fractions */
    double avg_alcroot;                 /* Average coarse root growing season allocation fractions */
    double avg_alleaf;                  /* Average leaf growing season allocation fractions */
    double avg_alroot;                  /* Average fine root growing season allocation fractions */
    double avg_alstem;                  /* Average stem growing season allocation fractions */
    double branch;                      /* branch c (t/ha) */
    double branchn;                     /* branch n (t/ha) */
    double canht;                       /* canopy height (m) */
    double croot;                       /* coarse root c (t/ha) */
    double crootn;                      /* coarse root n (t/ha) */
    double cstore;                      /* C store for deciduous model (t/ha) */
    double inorgn;                      /* Inorganic soil N pool - dynamic (t/ha) */
    /* phosphorus (pcycle) */
    double shootp; /* shoot P (t/ha) */
    double rootp; /* fine root P (t/ha) */
    double crootp; /* coarse root P (t/ha) */
    double branchp; /* branch P (t/ha) */
    double stemp; /* stem P = stempimm + stempmob (t/ha) */
    double stempimm; /* immobile stem P (t/ha) */
    double stempmob; /* mobile stem P (t/ha) */
    double pstore; /* P store, deciduous (t/ha) */
    double shootpc; /* shoot P:C */
    double rootpc; /* root P:C */
    double structsurfp; /* surface structural litter P (t/ha) */
    double structsoilp; /* soil structural litter P (t/ha) */
    double metabsurfp; /* surface metabolic litter P (t/ha) */
    double metabsoilp; /* soil metabolic litter P (t/ha) */
    double activesoilp; /* active SOM P (t/ha) */
    double slowsoilp; /* slow SOM P (t/ha) */
    double passivesoilp; /* passive SOM P (t/ha) */
    double inorglabp; /* labile inorganic P (t/ha) */
    double inorgsorbp; /* sorbed inorganic P (t/ha) */
    double inorgssorbp; /* strongly sorbed inorganic P (t/ha) */
    double inorgoccp; /* occluded P (t/ha) */
    double inorgparp; /* parent material P (t/ha) */
    double inorgavlp; /* available mineral P = labile + sorbed (t/ha) */
    double inorgp; /* total inorganic P (t/ha) */
    double plantp; /* plant P (t/ha) */
    double litterp; /* litter P (t/ha) */
    double litterpag; /* aboveground litter P (t/ha) */
    double litterpbg; /* belowground litter P (t/ha) */
    double soilp; /* soil organic + inorganic P (t/ha) */
    double totalp; /* total ecosystem P (t/ha) */
    double p_to_alloc_shoot; /* deciduous P allocation */
    double p_to_alloc_root;
    double p_to_alloc_croot;
    double p_to_alloc_branch;
    double p_to_alloc_stemimm;
    double p_to_alloc_stemmob;
    double lai;                         /* leaf area index m2 (leaf) m-2 (ground) */
    double lai_prescribed;   /* today's prescribed LAI (m2 m-2) */
    double fipar;
    double metabsoil;                   /* metabolic soil c (t/ha) */
    double metabsoiln;                  /* metabolic soil n (t/ha) */
    double metabsurf;                   /* metabolic surface c (t/ha) */
    double metabsurfn;                  /* metabolic surface n (t/ha) */
    double nstore;                      /* N store for deciduous model (t/ha) */
    double passivesoil;                 /* passive C som pool (t/ha) */
    double passivesoiln;                /* passive N som pool (t/ha) */
    double pawater_root;                /* plant available water - root zone (mm) */
    double pawater_topsoil;             /* plant available water - top soil(mm) */
    double prev_sma;
    double root;                        /* root c (t/ha) */
    double root_depth;                  /* rooting depth, Dmax (m) */
    double rootn;                       /* root n (t/ha) */
    double sapwood;
    double shoot;                       /* shoot c (t/ha) */
    double shootn;                      /* shoot n (t/ha) */
    double sla;                         /* specific leaf area */
    double slowsoil;                    /* slow C som pool (t/ha) */
    double slowsoiln;                   /* slow N som pool (t/ha) */
    double stem;
    double stemn;                       /* Stem N (t/ha) = stemnimm + stemnmob */
    double stemnimm;
    double stemnmob;
    double structsoil;                  /* soil structural c (t/ha) */
    double structsoiln;                 /* soil structural n (t/ha) */
    double structsurf;                  /* surface structural c (t/ha) */
    double structsurfn;                 /* surface structural n (t/ha) */
    double shootnc;
    double rootnc;
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
    double n_to_alloc_root;
    double c_to_alloc_root;
    double c_to_alloc_croot;
    double n_to_alloc_croot;
    double c_to_alloc_branch;
    double n_to_alloc_branch;
    double c_to_alloc_stem;
    double n_to_alloc_stemmob;
    double n_to_alloc_stemimm;
    double anpp;
    double litterc;
    double littern;
    double littercbg;
    double littercag;
    double litternag;
    double litternbg;
    double plantc;
    double plantn;
    double totaln;
    double totalc;
    double soilc;
    double soiln;
    double canopy_store;
    double psi_s_topsoil;
    double psi_s_root;

    /* hydraulics */
    double *thickness;
    double *root_mass;
    double *root_length;
    double *layer_depth;
    double *wetting_bot;
    double *wetting_top;
    double *water_frac;
    double initial_water;
    double weighted_swp;
    double dry_thick;   /* Thickness of dry soil layer above water table (m)*/
    int    rooted_layers;
    double predawn_swp;     /* MPa */
    double midday_lwp;     /* MPa */
    double lwp;
    double midday_plc;     /* midday loss of xylem conductance (%) */
} state;

typedef struct {
    double a0rhizo; /* minimum allocation to rhizodeposition [0.0-0.1] */
    double a1rhizo; /* slope of allocation to rhizodeposition [0.2-1] */
    double actncmax;                        /* Active pool (=1/3) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double actncmin;                        /* Active pool (=1/15) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double adapt;
    double ageold;                          /* Plant age when max leaf N C ratio is lowest */
    double ageyoung;                        /* Plant age when max leaf N C ratio is highest */
    double albedo;
    double alpha_c4;                        /* quantium efficiency for C4 plants has no Ci and temp dependancy, so if a fixed constant. */
    double alpha_j;                         /* initial slope of rate of electron transport, used in calculation of quantum yield. Value calculated by Belinda */
    double b_root;
    double b_topsoil;
    double bdecay;                          /* branch and large root turnover rate (1/yr) */
    double branch0;                         /* constant in branch-stem allometry (trees) */
    double branch1;                         /* exponent in branch-stem allometry */
    double bretrans;                        /* branch n retranslocation fraction */
    double c_alloc_bmax;                    /* allocation to branches at branch n_crit. If using allometric model this is the max alloc to branches */
    double c_alloc_bmin;                    /* allocation to branches at zero branch n/c. If using allometric model this is the min alloc to branches */
    double c_alloc_cmax;                    /* allocation to coarse roots at n_crit. If using allometric model this is the max alloc to coarse roots */
    double c_alloc_fmax;                    /* allocation to leaves at leaf n_crit. If using allometric model this is the max alloc to leaves */
    double c_alloc_fmin;                    /* allocation to leaves at zero leaf n/c. If using allometric model this is the min alloc to leaves */
    double c_alloc_rmax;                    /* allocation to roots at root n_crit. If using allometric model this is the max alloc to fine roots */
    double c_alloc_rmin;                    /* allocation to roots at zero root n/c. If using allometric model this is the min alloc to fine roots */
    double cfracts;                         /* conversion factor -> from dry matter (DM) to C */
    double crdecay;                         /* coarse roots turnover rate (1/yr) */
    double cretrans;                        /* coarse root n retranslocation fraction */
    double croot0;                          /* constant in coarse_root-stem allometry (trees) */
    double croot1;                          /* exponent in coarse_root-stem allometry */
    double ctheta_root;                     /* Fitted parameter based on Landsberg and Waring */
    double ctheta_topsoil;                  /* Fitted parameter based on Landsberg and Waring */
    double cue;                             /* carbon use efficiency, or the ratio of NPP to GPP */
    double d0;
    double d0x;                             /* Length scale for exponential decline of Umax(z) */
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
    double edv;                             /* Deactivation energy for Rubisco (J mol-1), <= 0: plain Arrhenius for Vcmax */
    double delsv;                           /* Entropy term for Rubisco deactivation (J mol-1 K-1) */
    double photo_tlow;                      /* Leaf/air T (deg C) below which Jmax/Vcmax = 0 */
    double photo_thigh;                     /* T (deg C) above which there is no low T reduction of Jmax/Vcmax */
    double faecescn;
    double faecesn;                         /* Faeces C:N ratio */
    double fdecay;                          /* foliage turnover rate (1/yr) */
    double fdecaydry;                       /* Foliage turnover rate - dry soil (1/yr) */
    double fhw;                             /* n:c ratio of stemwood - immobile pool and new ring */
    double finesoil;                        /* clay+silt fraction */
    double fix_lai;                       /* value to fix LAI to, control fixed_lai flag must be set */
    double nc_co2;                        /* netCDF forcing: CO2 (ppm) if the file has no CO2air */
    double nc_ndep;                       /* netCDF forcing: N deposition (t N ha-1 yr-1) */
    double nc_nfix;                       /* netCDF forcing: N fixation (t N ha-1 yr-1) */
    double fracfaeces;                      /* Fractn of grazd C that ends up in faeces (0..1) */
    double fracteaten;                      /* Fractn of leaf prodn eaten by grazers */
    double fractosoil;                      /* Fractn of grazed N recycled to soil:faeces+urine */
    double fractup_soil;                    /* fraction of uptake from top soil layer */
    double fretrans;                        /* foliage n retranslocation fraction - 46-57% in young E. globulus trees - see Corbeels et al 2005 ecological modelling 187, pg 463. Roughly 50% from review Aerts '96 */
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
    double kext;                            /* extinction coefficient */
    double kn;                              /* extinction coefficient of nitrogen in the canopy, assumed to be 0.3 by defaul which comes half from Belinda's head and is supported by fig 10 in Lloyd et al. Biogeosciences, 7, 1833–1859, 2010 */
    double knl;
    double ko25;                            /* Base rate for oxygenation by Rubisco at 25degC [umol mol-1]. Note value in Bernacchie 2001 is in mmol!! */
    double kq10;                            /* exponential coefficient for Rm vs T */
    double kr;                              /* N uptake coefficent (0.05 kg C m-2 to 0.5 tonnes/ha) see Silvia's PhD, Dewar and McM, 96. */
    double lad;                             /* Leaf angle distribution: 0 = spherical leaf angle distribution; 1 = horizontal leaves; -1 = vertical leaves */
    double lai_closed;                      /* LAI of closed canopy (max cover fraction is reached (m2 (leaf) m-2 (ground) ~ 2.5) */
    double latitude;                        /* latitude (degrees, negative for south) */
    double longitude;                       /* longitude (degrees, negative for west) */
    double leafsap0;                        /* leaf area  to sapwood cross sectional area ratio when Height = Height0 (mm^2/mm^2) */
    double leafsap1;                        /* leaf to sap area ratio when Height = Height1 (mm^2/mm^2) */
    double ligfaeces;                       /* Faeces lignin as fractn of biomass */
    double ligroot;                         /* lignin-to-biomass ratio in root litter; Values from White et al. = 0.22  - Value in Smith et al. 2013 = 0.16, note subtly difference in eqn C9. */
    double ligshoot;                        /* lignin-to-biomass ratio in leaf litter; Values from White et al. DBF = 0.18; ENF = 0.24l; GRASS = 0.09; Shrub = 0.15 - Value in smith et al. 2013 = 0.2, note subtly difference in eqn C9. */
    double liteffnc;
    double max_intercep_lai;                /* canopy LAI at which interception is maximised. */
    double measurement_temp;                /* temperature Vcmax/Jmax are measured at, typical 25.0 (celsius)  */
    double ncbnew;                          /* N alloc param: new branch N C at critical leaf N C */
    double ncbnewz;                         /* N alloc param: new branch N C at zero leaf N C */
    double nccnew;                          /* N alloc param: new coarse root N C at critical leaf N C */
    double nccnewz;                         /* N alloc param: new coarse root N C at zero leaf N C */
    double ncmaxfold;                       /* max N:C ratio of foliage in old stand, if the same as young=no effect */
    double ncmaxfyoung;                     /* max N:C ratio of foliage in young stand, if the same as old=no effect */
    double ncmaxr;                          /* max N:C ratio of roots */
    double ncrfac;                          /* N:C of fine root prodn / N:C c of leaf prodn */
    double ncwimm;                          /* N alloc param: Immobile stem N C at critical leaf N C */
    double ncwimmz;                         /* N alloc param: Immobile stem N C at zero leaf N C */
    double ncwnew;                          /* N alloc param: New stem ring N:C at critical leaf N:C (mob) */
    double ncwnewz;                         /* N alloc param: New stem ring N:C at zero leaf N:C (mobile) */
    double nf_crit;                         /* leaf N:C below which N availability limits productivity  */
    double nf_min;                          /* leaf N:C minimum N concentration which allows productivity */
    double nmax;
    double nmin;                            /* (bewdy) minimum leaf n for +ve p/s (g/m2) */
    double nmin0;                           /* mineral N pool corresponding to Actnc0,etc (g/m2) */
    /* phosphorus (pcycle) */
    double actpcmax; /* active SOM P:C max */
    double actpcmin; /* active SOM P:C min */
    double slowpcmax; /* slow SOM P:C max */
    double slowpcmin; /* slow SOM P:C min */
    double passpcmax; /* passive SOM P:C max */
    double passpcmin; /* passive SOM P:C min */
    double pmin0; /* labile P at which SOM P:C is at its min (g m-2) */
    double structcp; /* structural litter C:P */
    double structratp; /* structural litter P:C as a fraction of metabolic (strpfloat) */
    double pcmaxfyoung; /* max leaf P:C, young stand */
    double pcmaxfold; /* max leaf P:C, old stand */
    double pcrfac; /* fine root P:C as a fraction of leaf */
    double pcbnew; /* new branch P:C, old */
    double pcbnewz; /* new branch P:C, young */
    double pccnew; /* new coarse root P:C, old */
    double pccnewz; /* new coarse root P:C, young */
    double pcwimm; /* immobile stem P:C, old */
    double pcwimmz; /* immobile stem P:C, young */
    double pcwnew; /* new stem P:C, old */
    double pcwnewz; /* new stem P:C, young */
    double pf_min; /* min leaf P:C */
    double fretransp; /* leaf P retranslocation fraction */
    double prescribed_leaf_pc; /* leaf P:C when the P cycle is off (photosynthesis) */
    double faecescp; /* faeces C:P */
    double fractosoilp; /* fraction of grazed P returned to soil */
    double prateuptake; /* P uptake rate constant (yr-1) */
    double krp; /* P uptake root coefficient (t C/ha) */
    double puptakez; /* constant P uptake (t P/ha/yr, puptake_model 0) */
    double prateloss; /* labile P leaching rate (yr-1) */
    double smax; /* max sorbed P (t P/ha; 700 g m-2 Yang et al. 2016) */
    double ks; /* Langmuir half-saturation of sorption (t P/ha; 0.5 g m-2) */
    double rate_sorb_ssorb; /* sorbed -> strongly sorbed P (yr-1) */
    double rate_ssorb_occ; /* strongly sorbed -> occluded P (yr-1) */
    double psecmnp; /* strongly sorbed -> mineral P, constant (d-1; text_effect_p 0) */
    double phmin; /* CENTURY pH range, min */
    double phmax; /* CENTURY pH range, max */
    double phtextmin; /* CENTURY ssorb -> mineral rate, min (d-1) */
    double phtextmax; /* CENTURY ssorb -> mineral rate, max (d-1) */
    double phtextslope; /* CENTURY sand slope */
    double soilph; /* soil pH */
    double p_atm_deposition; /* atmospheric P deposition (t P/ha/yr) */
    double p_rate_par_weather; /* parent material P weathering rate (yr-1) */
    double max_p_biochemical; /* max biochemical P mineralisation (t P/ha/yr; Wang et al. 2007) */
    double crit_n_cost_of_p; /* N cost of P uptake above which phosphatase is made (g N/g P) */
    double biochemical_p_constant; /* biochemical P mineralisation half-saturation (g N/g P) */
    double kp_canopy; /* extinction coefficient of leaf P in the canopy */
    double jmaxpa; /* Jmax vs leaf P (Walker et al. 2014) */
    double jmaxpb;
    double vcmaxpa; /* Vcmax vs leaf P (Walker et al. 2014) */
    double vcmaxpb;
    double nmincrit;                        /* Critical mineral N pool at max soil N:C (g/m2) (Parton et al 1993, McMurtrie et al 2001). */
    double ntheta_root;                     /* Fitted parameter based on Landsberg and Waring */
    double ntheta_topsoil;                  /* Fitted parameter based on Landsberg and Waring */
    double nuptakez;                        /* constant N uptake per year (1/yr) */
    double oi;                              /* intercellular concentration of O2 [umol mol-1] */
    double passivesoilnz;
    double passivesoilz;
    double passncmax;                       /* Passive pool (=1/7) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double passncmin;                       /* Passive pool (=1/10) N:C of new SOM - when Nmin=Nmin0 [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double prescribed_leaf_NC;              /* If the N-Cycle is switched off this needs to be set, e.g. 0.03 */
    double previous_ncd;                    /* In the first year we don't have last years data, so I have precalculated the average of all the november-jan chilling values  */
    double psi_sat_root;                    /* MPa */
    double psi_sat_topsoil;                 /* MPa */
    double psie_topsoil;                    /* Soil water potential at saturation (m) */
    double psie_root;                       /* Soil water potential at saturation (m) */
    double qs;                              /* exponent in water stress modifier, =1.0 JULES type representation, the smaller the values the more curved the depletion.  */
    double r0;                              /* root C at half-maximum N uptake (kg C/m3) */
    double rateloss;                        /* Rate of N loss from mineral N pool (/yr) */
    double rateuptake;                      /* ate of N uptake from mineral N pool (/yr) from here? http://face.ornl.gov/Finzi-PNAS.pdf */
    double rdecay;                          /* root turnover rate (1/yr) */
    double rdecaydry;                       /* root turnover rate - dry soil (1/yr) */
    double resp_coeff;                      /* Respiration rate: from LPJ ENF, EBF, C3G = 1.2, Trop EBF, C4G = 0.2 */
    double retransmob;                      /* Fraction stem mobile N retranscd (/yr) */
    double rfmult;
    double rooting_depth;                   /* Rooting depth (mm) */
    char   rootsoil_type[STRING_LENGTH];
    double rretrans;                        /* root n retranslocation fraction */
    double sapturnover;                     /* Sapwood turnover rate: conversion of sapwood to heartwood (1/yr) */
    double sla;                             /* specific leaf area (m2 one-sided/kg DW) */
    double slamax;                          /* (if equal slazero=no effect) specific leaf area new fol at max leaf N/C (m2 one-sided/kg DW) */
    double slazero;                         /* (if equal slamax=no effect) specific leaf area new fol at zero leaf N/C (m2 one-sided/kg DW) */
    double slowncmax;                       /* Slow pool (=1/15) N:C ratio of new SOM - maximum [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double slowncmin;                       /* Slow pool (=1/40) N:C of new SOM - when Nmin=Nmin0" [units: gN/gC]. Based on forest version of CENTURY (Parton et al. 1993), see Appendix, McMurtrie 2001, Tree Physiology. */
    double store_transfer_len;
    double structcn;                        /* C:N ratio of structural bit of litter input */
    double structrat;                       /* structural input n:c as fraction of metab */
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
    double fmleaf;
    double fmroot;
    double decayrate[7];
    double fmfaeces;
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

    /* hydraulics */
    double layer_thickness;                 /* Soil layer thickness (m) */
    int    soil_layers;                     // Number of soil layers
    int    core;                            // Layer below soil layers, i.e. 20 + 1
    double root_k;    /* mass of roots for reaching 50% maximum depth (g m-2) */
    double root_radius;  /* (m) */
    double root_density; /* g biomass m-3*/
    double max_depth;    /* (m) */
    double root_resist;
    double gs_min;

    /* gs_opt (Sperry profit maximisation), hydraulics */
    double root_psi_crit;    /* water potential at which roots stop taking
                                up water, sets the layer uptake weights
                                (MPa), JULES root_psi_crit (was min_lwp) */
    double kp;               /* maximum plant (xylem) hydraulic conductance
                                per unit leaf area (mmol m-2 s-1 MPa-1),
                                JULES kmax_pft */
    double p50;              /* xylem water potential at 50% loss of
                                conductance (MPa) */
    double p88;              /* ... and at 88% (MPa) */
    double kcrit_frac;       /* loss of conductance at which the xylem fails,
                                kcrit = (1 - kcrit_frac) kp (-) */
    double gs_opt_gl_max;    /* max leaf conductance for H2O (m s-1), <= 0
                                for no cap, JULES som_gl_max */
    int    gs_opt_n_sample;  /* Ci samples, flat search */
    int    gs_opt_n_prescan; /* Ci samples, golden search prescan */
    int    gs_opt_n_golden;  /* golden section iterations */
    double gs_opt_peak_e;    /* daily gs_opt: the hydraulic cost is on the
                                peak (midday) E = this x the am/pm mean E */
    /* plant_segments = segmented: share of the whole plant resistance and
       vulnerability of the root, stem & leaf (P50/P88 < -900 = the whole
       plant values) */
    double seg_frac_root, seg_frac_stem, seg_frac_leaf;
    double p50_root, p50_stem, p50_leaf;
    double p88_root, p88_stem, p88_leaf;

    /* not shared via cmd line */
    double *potA;
    double *potB;
    double *cond1;
    double *cond2;
    double *cond3;
    double *porosity;
    double *field_capacity;
    /* van Genuchten / Brooks-Corey soil (JULES conventions & names; used for
       all layers): VG n = 1 + 1/b, alpha = 1/sathh, Mualem L = 0.5 */
    double soil_b;           /* Brooks-Corey b, or 1/(n-1) for VG (-) */
    double soil_sathh;       /* |psi| at saturation (BC) or 1/alpha (VG) (m) */
    double soil_satcon;      /* saturated conductivity (m s-1) */
    double soil_sm_sat;      /* saturated water content (m3 m-3) */
    double soil_sm_res;      /* residual water content, VG (m3 m-3) */
    double soil_psi_open;    /* bound_soil_psi: max soil water potential (MPa) */
    double soil_psi_close;   /* bound_soil_psi: min / dry soil matching point (MPa) */
    double ds_psi;           /* dry soil correction: psi at zero water content (MPa) */
    double ds_min_depth;     /* dry soil correction only below this depth (m) */
    /* soil evaporation: JULES (gs_nvg, gsoil_f) & Or/Decker (CABLE psm) */
    double gs_nvg;           /* JULES bare soil conductance (m s-1) */
    double gsoil_f;          /* JULES scaling of gsoil under the canopy (-) */
    double soil_tortuosity;  /* SPA (gday) soil evaporation: tortuosity of
                                the dry surface layer (-) */
    double or_sublayer_dz;   /* Or: viscous sublayer thickness (m) */
    double litter_dz_per_c;  /* Or: litter depth per litter C (m per t C ha-1) */
    double litter_c;         /* Or: surface litter C (t C ha-1), <0 = simulated */

    /* two-leaf canopy radiation & leaf boundary layer (CABLE) */
    double leaf_tau_vis;     /* leaf transmittance, visible (-) */
    double leaf_tau_nir;     /* leaf transmittance, near infrared (-) */
    double leaf_refl_vis;    /* leaf reflectance, visible (-) */
    double leaf_refl_nir;    /* leaf reflectance, near infrared (-) */
    double leaf_chi;         /* leaf angle distribution, 0 = spherical (-) */
    double soil_refl;        /* soil reflectance, SW mean (-) */
    double shelter;          /* leaf boundary layer sheltering factor (-) */
    double wind_height;      /* height of the wind forcing (m), <= canht
                                means the wind is at the canopy top */

    /* JULES canopy interception (canopy_evap_model = jules) */
    double catch0;           /* canopy capacity with no leaves (mm) */
    double dcatch_dlai;      /* canopy capacity per unit LAI (mm) */
    double rain_area_frac;   /* fraction of the area rain falls on (-) */
    int     wetting;         /* number of wetting layers */


} params;

typedef struct {

    double *year;
    double *rain;
    double *par;
    double *tair;
    double *tsoil;
    double *co2;
    double *ndep;
    double *nfix;       /* N inputs from biological fixation (t/ha/timestep (d/30min)) */
    double *lai;        /* prescribed LAI (m2 m-2), only if prescribed_lai */
    double *lwdown;     /* downward LW (W m-2), NULL if not in the forcing */
    double *lwdown_am, *lwdown_pm;  /* daily: daylight am/pm means */
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
    double nfix;       /* N inputs from biological fixation (t/ha/timestep (d/30min)) */
    double tsoil;
    double lwdown;     /* downward LW (W m-2), < 0 if not in the forcing */
    double lwdown_am, lwdown_pm;  /* daily am/pm means, < 0 if missing */

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
    double gsc_am, gsc_pm;  /* daily gs_opt: canopy gs for CO2 chosen by
                               MATE (mol m-2 s-1), used by the water
                               balance */
    double gpp_pm;
    double gpp;
    double npp;
    double nep;
    double auto_resp;
    double hetero_resp;
    double retrans;
    double apar;

    /* n */
    double nuptake;
    double nloss;
    /* phosphorus (pcycle, t P/ha/d) */
    double retransp;
    double leafretransp;
    double puptake;
    double ploss;
    double pgross;
    double pimmob;
    double plittrelease;
    double pmineralisation;
    double ppleaf;
    double pproot;
    double ppcroot;
    double ppbranch;
    double ppstemimm;
    double ppstemmob;
    double deadleafp;
    double deadrootp;
    double deadcrootp;
    double deadbranchp;
    double deadstemp;
    double peaten;
    double purine;
    double p_surf_struct_litter;
    double p_surf_metab_litter;
    double p_soil_struct_litter;
    double p_soil_metab_litter;
    double p_surf_struct_to_slow;
    double p_soil_struct_to_slow;
    double p_surf_struct_to_active;
    double p_soil_struct_to_active;
    double p_surf_metab_to_active;
    double p_soil_metab_to_active;
    double p_active_to_slow;
    double p_active_to_passive;
    double p_slow_to_active;
    double p_slow_to_passive;
    double p_slow_biochemical;
    double p_passive_to_active;
    double p_lab_in;
    double p_lab_out;
    double p_sorb_in;
    double p_sorb_out;
    double p_min_to_ssorb;
    double p_ssorb_to_min;
    double p_ssorb_to_occ;
    double p_par_to_min;
    double p_atm_dep;
    double lprate;
    double rprate;
    double bprate;
    double cprate;
    double wpimrate;
    double wpmobrate;
    double npassive;        /* n passive -> active */
    double ngross;          /* N gross mineralisation */
    double nimmob;          /* N immobilisation in SOM */
    double nlittrelease;    /* N rel litter = struct + metab */
    double activelossf;     /* frac of active C -> CO2 */
    double nmineralisation;

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
    double day_ppt;
    double day_wbal;

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
    double deadsapwood;

    /* grazing stuff */
    double ceaten;          /* C consumed by grazers (t C/ha/y) */
    double neaten;          /* N consumed by grazers (t C/ha/y) */
    double faecesc;         /* Flux determined by faeces C:N */
    double nurine;          /* Rate of N input to soil in urine (t/ha/y) */

    double leafretransn;

    /* C&N Surface litter */
    double surf_struct_litter;
    double surf_metab_litter;
    double n_surf_struct_litter;
    double n_surf_metab_litter;

    /* C&N Root Litter */
    double soil_struct_litter;
    double soil_metab_litter;
    double n_soil_struct_litter;
    double n_soil_metab_litter;

    /* C&N litter fluxes to slow pool */
    double surf_struct_to_slow;
    double soil_struct_to_slow;
    double n_surf_struct_to_slow;
    double n_soil_struct_to_slow;

    /* C&N litter fluxes to active pool */
    double surf_struct_to_active;
    double soil_struct_to_active;
    double n_surf_struct_to_active;
    double n_soil_struct_to_active;

    /* Metabolic fluxes to active pool */
    double surf_metab_to_active;
    double soil_metab_to_active;
    double n_surf_metab_to_active;
    double n_soil_metab_to_active;

    /* C fluxes out of active pool */
    double active_to_slow;
    double active_to_passive;
    double n_active_to_slow;
    double n_active_to_passive;

    /* C&N fluxes from slow to active pool */
    double slow_to_active;
    double slow_to_passive;
    double n_slow_to_active;
    double n_slow_to_passive;

    /* C&N fluxes from passive to active pool */
    double passive_to_active;
    double n_passive_to_active;

    /* C & N source fluxes from the active, slow and passive pools */
    double c_into_active;
    double c_into_slow;
    double c_into_passive;

    /* CO2 flows to the air */
    double co2_to_air[7];

    /* C allocated fracs  */
    double alleaf;
    double alroot;
    double alcroot;
    double albranch;
    double alstem;

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

    double lnrate;
    double bnrate;
    double wnimrate;
    double wnmobrate;
    double cnrate;

    /* priming/exudation */
    double root_exc;
    double root_exn;
    double co2_released_exud;
    double factive;
    double rtslow;
    double rexc_cue;

    double ninflow;

    /* hydraulics */
    double *soil_conduct;
    double *swp;
    double *soilR;
    double *fraction_uptake;
    double *ppt_gain;
    double *water_loss;
    double *water_gain;
    double *est_evap;
    double total_soil_resist;

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
    double lwp_leaf[2];     /* leaf water potential (MPa) */
    double fwsoil_leaf[2];  /* gs_opt water stress factor, A / A(wet soil) */
    double an_canopy;       /* canopy net photosynthesis (umol m-2 s-1) */
    double rd_canopy;       /* canopy respiration in the light (umol m-2 s-1) */
    double gsc_canopy;      /* canopy stomatal conductance to CO2 (mol m-2 s-1) */
    double apar_canopy;     /* canopy abs photosyn. active rad. (umol m-2 s-1) */
    double omega_canopy;    /* canopy decoupling coefficient (-) */
    double trans_canopy;    /* canopy transpiration (mm 30min-1) */
    double rnet_canopy;     /* canopy net radiation (W m-2) */
    double lwp_canopy;      /* Leaf water potential for the canopy(MPa) */
    double N0;              /* top of canopy nitrogen (g N m-2)) */
    double elevation;       /* sun elevation angle in degrees */
    double cos_zenith;      /* cos(zenith angle of sun) in radians */
    double diffuse_frac;    /* Fraction of incident rad which is diffuse (-) */
    double direct_frac;     /* Fraction of incident rad which is beam (-) */
    double tleaf_new;       /* new leaf temperature (deg C) */
    double dleaf;           /* leaf VPD (Pa) */
    double Cs;              /* CO2 conc at the leaf surface (umol mol-1) */
    double kb;              /* beam radiation ext coeff of canopy */
    double kd;              /* diffuse radiation ext coeff of canopy */
    double gradis[2];       /* big-leaf radiative conductance (mol m-2 s-1) */
    double gbhu[2];         /* big-leaf forced convection boundary layer
                               conductance for heat (mol m-2 s-1) */
    double qssabs;          /* SW absorbed by the soil (W m-2) */
    double rnet_soil;       /* soil net radiation (W m-2) */
    double scalex[2];      /* scale from single leaf to canopy */
    double *cz_store;       /* Array to hold coz zenith angles */
    double *ele_store;      /* Array to hold elevations */
    double *df_store;       /* Array to hold diffuse fractions */

    /* gs_opt */
    double kl_leaf[2];      /* xylem conductance at the leaf water potential
                               (mmol m-2 s-1 MPa-1, per unit leaf area) */
    double kl_canopy;       /* ... canopy mean */
    double psi_stem_leaf[2]; /* stem water potential feeding each big leaf
                               (MPa, = leaf when not segmented) */
    double psi_stem_canopy;
    double e_supply;        /* transpiration the rooted soil can supply
                               this step, water above theta(root_psi_crit)
                               (mmol m-2 s-1, ground), < 0 no limit */
    double tair_canopy;     /* canopy air temperature (deg C) and VPD (Pa)
                               the leaves saw */
    double vpd_canopy;
    double dtc_prev, dec_prev;  /* canopy air - reference T (K) and vapour
                                   pressure (Pa) of the last daytime step,
                                   to start the canopy air iteration */
    double cs_leaf[2], dleaf_leaf[2];  /* each big leaf's Cs and VPD, kept
                                          between canopy air passes */
    double gbv_leaf[2];     /* big leaf boundary layer conductance for H2O
                               from the last energy balance (mol m-2 s-1) */

} canopy_wk;


typedef struct {
    double  *ystart;
    double   *yscal;
    double   *y;
    double   *dydx;
    double   *xp;
	double  **yp;
    int       N;
    int       kmax;

    double   *ak2;
    double   *ak3;
    double   *ak4;
    double   *ak5;
    double   *ak6;
    double   *ytemp;
    double   *yerr;



} nrutil;

/* running sums over the last spin-up cycle, for the SAS spin-up */
typedef struct {
    long   ndays;
    double dr[7];                   /* litter & SOM decay rates (d-1) */
    double surf_struct_litter;      /* litter C inputs (t C ha-1 d-1) */
    double surf_metab_litter;
    double soil_struct_litter;
    double soil_metab_litter;
    double cpbranch;                /* woody growth & turnover */
    double cpstem;
    double cpcroot;
    double deadbranch;
    double deadstems;
    double deadcroots;
    double metablsoil_nc;           /* pool N:C ratios */
    double metabsurf_nc;
    double structsoil_nc;
    double structsurf_nc;
    double activesoil_nc;
    double slowsoil_nc;
    double passivesoil_nc;
} fast_spinup;

#endif
