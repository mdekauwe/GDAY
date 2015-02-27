#ifndef CONSTANTS_H
#define CONSTANTS_H
/*
    constants.h
        - Holds all the macro definitions for constants used in GDAY
*/


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
#define MM_TO_M  0.001
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
#define KG_AS_TONNES 1E-3
#define TONNES_AS_KG 1.0 / KG_AS_TONNES



#endif /* CONSTANTS */
