#ifndef GDAY_H
#define GDAY_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>


#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define STRING_LENGTH 2000





#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CLIP(x) ((x)<0. ? 0. : ((x)>1. ? 1. : (x)))

/* Stomatal conductanct models */
#define MEDLYN 0

/* Photosynthesis models */
#define MATE 0
#define BEWDY 1

/* Respiration models */
#define FIXED 0
#define TEMPERATURE 1
#define BIOMASS 2

/* Allocation models */
#define FIXED 0
#define GRASSES 1
#define ALLOMETRIC 2

/* Ps pathway */
#define C3 0
#define C4 1

#define DAILY 0
#define END 1

typedef struct {
    FILE *ifp;
    FILE *ofp;
    char  cfg_fname[STRING_LENGTH];
    char  met_fname[STRING_LENGTH];
    char  out_fname[STRING_LENGTH];
    char  out_param_fname[STRING_LENGTH];
    char  git_hash[STRING_LENGTH];
    int   adjust_rtslow;
    int   alloc_model;
    int   assim_model;
    int   calc_sw_params;
    int   deciduous_model;
    int   disturbance;
    int   exudation;
    int   fixed_stem_nc;
    int   fixleafnc;
    int   grazing;
    int   gs_model;
    int   hurricane;
    int   model_optroot;
    int   modeljm;
    int   ncycle;
    int   num_years;
    int   nuptake_model;
    int   output_ascii;
    int   passiveconst;
    int   print_options;
    int   ps_pathway;
    int   respiration_model;
    int   strfloat;
    int   sw_stress_model;
    int   trans_model;
    int   use_eff_nc;
    int   water_stress;
    int   num_days;
    char  git_code_ver[STRING_LENGTH];
    int   spin_up;
    int   PRINT_GIT;
} control;


typedef struct {
    double activesoil;
    double activesoiln;
    double age;
    double avg_albranch;
    double avg_alcroot;
    double avg_alleaf;
    double avg_alroot;
    double avg_alstem;
    double branch;
    double branchn;
    double canht;
    double croot;
    double crootn;
    double cstore;
    double inorgn;
    double lai;
    double fipar;
    double max_lai;
    double max_shoot;
    double metabsoil;
    double metabsoiln;
    double metabsurf;
    double metabsurfn;
    double nstore;
    double passivesoil;
    double passivesoiln;
    double pawater_root;
    double pawater_topsoil;
    double prev_sma;
    double root;
    double root_depth;
    double rootn;
    double sapwood;
    double shoot;
    double shootn;
    double sla;
    double slowsoil;
    double slowsoiln;
    double stem;
    double stemn;
    double stemnimm;
    double stemnmob;
    double structsoil;
    double structsoiln;
    double structsurf;
    double structsurfn;
    double shootnc;
    double rootnc;
    double remaining_days[366];
    double wtfac_root;
    double wtfac_topsoil;
    double delta_sw_store;
    double grw_seas_stress;
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


} state;

typedef struct {
    double a1;
    double ac;
    double actncmax;
    double actncmin;
    double adapt;
    double ageold;
    double ageyoung;
    double albedo;
    double alpha_c4;
    double alpha_j;
    double b_root;
    double b_topsoil;
    double bdecay;
    double branch0;
    double branch1;
    double bretrans;
    double burn_specific_yr;
    double c_alloc_bmax;
    double c_alloc_bmin;
    double c_alloc_cmax;
    double c_alloc_fmax;
    double c_alloc_fmin;
    double c_alloc_rmax;
    double c_alloc_rmin;
    double cfracts;
    double crdecay;
    double cretrans;
    double croot0;
    double croot1;
    double ctheta_root;
    double ctheta_topsoil;
    double cue;
    double d0;
    double d0x;
    double d1;
    double delsj;
    double density;
    double direct_frac;
    double displace_ratio;
    double disturbance_doy;
    double dz0v_dh;
    double eac;
    double eag;
    double eaj;
    double eao;
    double eav;
    double edj;
    double faecescn;
    double faecesn;
    double fdecay;
    double fdecaydry;
    double fhw;
    double finesoil;
    double fracfaeces;
    double fracteaten;
    double fractosoil;
    double fractup_soil;
    double fretrans;
    double g1;
    double gamstar25;
    double growth_efficiency;
    double height0;
    double height1;
    double heighto;
    double htpower;
    double hurricane_doy;
    double hurricane_yr;
    double intercep_frac;
    double jmax;
    double jmaxna;
    double jmaxnb;
    double jv_intercept;
    double jv_slope;
    double kc25;
    double kdec1;
    double kdec2;
    double kdec3;
    double kdec4;
    double kdec5;
    double kdec6;
    double kdec7;
    double kext;
    double knl;
    double ko25;
    double kq10;
    double kr;
    double lai_closed;
    double latitude;
    double leafsap0;
    double leafsap1;
    double ligfaeces;
    double ligroot;
    double ligshoot;
    double liteffnc;
    double max_intercep_lai;
    double measurement_temp;
    double ncbnew;
    double ncbnewz;
    double nccnew;
    double nccnewz;
    double ncmaxfold;
    double ncmaxfyoung;
    double ncmaxr;
    double ncrfac;
    double ncwimm;
    double ncwimmz;
    double ncwnew;
    double ncwnewz;
    double nf_crit;
    double nf_min;
    double nmax;
    double nmin;
    double nmin0;
    double nmincrit;
    double ntheta_root;
    double ntheta_topsoil;
    double nuptakez;
    double oi;
    double passivesoilnz;
    double passivesoilz;
    double passncmax;
    double passncmin;
    double prescribed_leaf_NC;
    double previous_ncd;
    double prime_y;
    double prime_z;
    double psi_sat_root;
    double psi_sat_topsoil;
    double qs;
    double r0;
    double rateloss;
    double rateuptake;
    double rdecay;
    double rdecaydry;
    double retransmob;
    double return_interval;
    double rfmult;
    double root_exu_CUE;
    double rooting_depth;
    char   rootsoil_type[STRING_LENGTH];
    double rretrans;
    double sapturnover;
    double sla;
    double slamax;
    double slazero;
    double slowncmax;
    double slowncmin;
    double store_transfer_len;
    double structcn;
    double structrat;
    double targ_sens;
    double theta;
    double theta_sat_root;
    double theta_sat_topsoil;
    double topsoil_depth;
    char   topsoil_type[STRING_LENGTH];
    double vcmax;
    double vcmaxna;
    double vcmaxnb;
    double watdecaydry;
    double watdecaywet;
    double wcapac_root;
    double wcapac_topsoil;
    double wdecay;
    double wetloss;
    double wretrans;
    double z0h_z0m;
    double fmleaf;
    double fmroot;
    double decayrate[7];
    double fmfaeces;
    int growing_seas_len;
} params;

typedef struct {
    double *year;
    double *prjday;
    double *sw_rad;
    double *tair;
    double *rain;
    double *tsoil;
    double *tam;
    double *tpm;
    double *vpd_am;
    double *vpd_pm;
    double *vpd_avg;
    double *co2;
    double *ndep;
    double *wind;
    double *atmos_press;
    double *par;
    double *wind_am;
    double *wind_pm;
    double *sw_rad_am;
    double *sw_rad_pm;
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
    double retrans;
    double apar;

    /* n */
    double nuptake;
    double nloss;
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
    double erain;
    double interception;
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
    double microbial_resp;
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
} fluxes;





void    read_met_data(char **, control *, met *);
int     calc_days_in_year(int);

int parse_ini_file(control *, params *, state *);
void usage(char **);


void run_sim(control *, fluxes *, met *, params *p, state *);
void spin_up_pools(control *, fluxes *, met *, params *p, state *);

/* gday main */
void   correct_rate_constants(params *, int output);
void   clparser(int, char **, control *);

void   reset_all_n_pools_and_fluxes(fluxes *, state *);
void   zero_stuff(control *, state *);
void   day_end_calculations(control *, params *, state *, int, int);

#endif /* GDAY_H */
