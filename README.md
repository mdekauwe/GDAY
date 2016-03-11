# GDAY model

GDAY (Generic Decomposition And Yield) is a simple ecosystem model that represents carbon, nitrogen, and water dynamics at the stand scale. The model can be run at either a daily time step, or sub-daily (i.e. 30-minutes). When the model is run at the sub-daily timescale, photosynthesis is calculated using a two-leaf approximation, otherwise photosynthesis is calculated following Sands (1995,1996). The sub-daily approach (photosynthesis & leaf energy balance) mirrors [MAESTRA](http://maespa.github.io/manual.html), without the complexity of the radiation treatment.

## Dependancies
The model is coded entirely in C without any dependancies. The wrapper files
for the example scripts and the script used to change GDAY parameter options,
are written in python. The old python version is still [online](https://github.com/mdekauwe/pygday).

## Installation
There is a Makefile in the src directory...

```bash
$ make clean ; make
```

I haven't done anything about an installation directory so you will need to move
the executable yourself.

## Running the model
A simple model usage can be displayed by calling GDAY as follows ::

```bash
$ gday -u
```

Alternative to run GDAY:

```bash
$ gday -p param_file.cfg
```

## Key References
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. *Ecological Applications*, 3, 666-681.
2. Medlyn, B. E., McMurtrie, R. E., Dewar, R. C. and Jeffreys, M. P. (2000), Soil processes dominate the long-term response of forest net primary productivity to increased temperature and atmospheric CO2 concentration, *Canadian Journal of Forest Research*, 30, 873–888.
3. Sands PJ (1995) Modelling canopy production. II. From single-leaf photosynthetic parameters to daily canopy photosynthesis. Australian Journal of Plant Physiology 22, 603-614.
4. Sands PJ (1996) Modelling canopy production. III. Canopy light-utilisation efficiency and its sensitivity to physiological environmental variables. Australian Journal of Plant Physiology 23, 103-114.

## Contacts
* [Martin De Kauwe](http://mdekauwe.github.io/).
* [Belinda Medlyn](<http://bio.mq.edu.au/people/person.php?user=bmedlyn).
