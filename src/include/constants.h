#ifndef CONSTANTS_H
#define CONSTANTS_H
/*
    constants.h
        - Holds all the macro definitions for constants used in GDAY
*/

#define MOL_2_MMOL 1000.0
#define MMOL_2_MOL 1E-03
#define NDAYS_IN_YR 365.25
#define DEG_TO_KELVIN 273.15
#define RGAS 8.314
#define SECS_IN_HOUR 3600.0
#define UMOL_TO_MOL 1E-6
#define MOL_C_TO_GRAMS_C 12.0
#define G_AS_TONNES 1E-6
#define M2_AS_HA 1E-4
#define KG_AS_G 1E+3
#define WATT_HR_TO_MJ 0.0036
#define MJ_TO_WATT_HR 1.0 / 0.0036
#define MM_TO_M  0.001
#define M_TO_MM  1000.0
#define GRAMS_C_TO_MOL_C 1.0 / 12.0
#define UMOL_TO_MOL 1E-6
#define MOL_TO_UMOL 1E6
#define GRAM_C_2_TONNES_HA G_AS_TONNES / M2_AS_HA
#define UMOL_2_GRAMS_C  UMOL_TO_MOL * MOL_C_TO_GRAMS_C
#define TONNES_HA_2_KG_M2 0.1
#define TONNES_HA_2_G_M2 100.0
#define YRS_IN_DAYS 1.0 / 365.25
#define DAYS_IN_YRS 365.25
#define G_M2_2_TONNES_HA 0.01
#define KG_M2_2_TONNES_HA 10.0
#define KPA_2_MPA 0.001
#define KG_AS_TONNES 1E-3
#define G_TO_KG 0.001
#define TONNES_AS_KG 1.0 / KG_AS_TONNES
#define CP 1010.0              /* specific heat of dry air (j kg-1 k-1) */
#define MASS_AIR 29.0E-3       /* mol mass air (kg mol-1) */
#define H2OLV0 2.501E6         /* latent heat H2O (J kg-1) */
#define H2OMW 18E-3            /* mol mass H20 (kg mol-1) */
#define DHEAT 21.5E-6          /* molecular diffusivity for heat */
#define GBVGBH 1.075           /* Ratio of Gbw:Gbh */
#define GSVGSC 1.57            /* Ratio of Gsw:Gsc */
#define GBHGBC 1.32            /* Ratio of Gbh:Gbc */
#define SIGMA 5.67e-8          /* Steffan Boltzman constant (W/m2/K4) */
#define LEAF_EMISSIVITY 0.95   /* Emissivity of thermal radiation by leaf */
#define KPA_2_PA 1000.
#define KPA_2_MPA 0.001
#define METER_OF_HEAD_TO_MPA 9.81 * KPA_2_MPA /* Height (m) x gravity (m/s2) = pressure (kPa) */
#define PA_2_KPA 0.001
#define CM_2_M 0.01

/* Solar radiaiton 1 W m-2 ~ 2.3 umol m-2 s-1 PAR
   Landsberg and Sands, Cp2, pg 20. (1.0 / 2.3) */
#define SW_2_PAR 2.3
#define PAR_2_SW 1.0 / SW_2_PAR
#define J_TO_MJ  1.0E-6
#define MJ_TO_J  1.0 / J_TO_MJ
#define J_2_UMOL 4.57               /* Conversion from J to umol quanta */
#define UMOL_2_JOL 1.0 / J_2_UMOL   /* Conversion from umol quanta to J */
#define SEC_2_HLFHR 1800.;
#define MOLE_WATER_2_G_WATER 18.02 /* oxygen = 16g/mol, hydrogren = 1.01 g/mol*/
#define SUNLIT 0
#define SHADED 1
#define NUM_LEAVES 2

#endif /* CONSTANTS */
