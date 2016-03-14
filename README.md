# GDAY model

GDAY (Generic Decomposition And Yield) is a simple ecosystem model that simulates carbon, nitrogen, and water dynamics at the stand scale. The model can be run at either a daily time step, or sub-daily (i.e. 30-minutes). When the model is run at the sub-daily timescale, photosynthesis is calculated using a two-leaf approximation (de Pury and Farquhar, 1997; Wang and Leuning, 1998), otherwise photosynthesis is calculated following Sands (1995;1996). The sub-daily approach (photosynthesis & leaf energy balance) mirrors [MAESTRA](http://maespa.github.io/manual.html), without the complexity of the radiation treatment.

**NOTE** the sub-daily version of the model is still in the development branch, but will be pulled into the main code asap.

<p style="text-align:center"><img src="doc/outline.png" width="500"/></p>

## Installation
The model is coded entirely in C without any dependancies. The wrapper files
for the example scripts and the script used to change parameter options,
are written in python. The old python version is still [online](https://github.com/mdekauwe/pygday).

There is a Makefile in the src directory...

```bash
$ make clean ; make
```

The Makefile will need to be edited by hand to set the $ARCH flag, which sets the installation path. Currently it is hardwired to my computer.

## Running the model
A simple model usage can be displayed by calling GDAY as follows:

```bash
$ gday -u
```

Presently, there are only two options which the user can set via the command line. All model options, including *all* of the model parameters are customisable via the parameter file.

To spin-up GDAY:

```bash
$ gday -s -p param_file.cfg
```

To run GDAY:

```bash
$ gday -p param_file.cfg
```

When the model is run it expects to find its "model state" (i.e. from a previous spin-up) in the parameter file. This state is automatically written the parameter file after the initial spin-up when the "print_options" flag has been set to "end", rather than "daily".


## Parameter file

GDAY expects a parameter file to be supplied as an argument (-p filename) on the command line. Parameter files follow the standard [.ini](https://en.wikipedia.org/wiki/INI_file) format, although only a relatively simple implementation has been coded into GDAY.

Parameter files are broken down into 6 section, namely [git], [files], [params], [control], [state] and [print]. The order of these sections shouldn't make any difference. The basic element contained in the parameter file is the key or property. Every key has a name and a value, delimited by an equals sign (=). The name appears to the left of the equals sign.

```ini
eac = 79430.0
```

As mentioned above, keys are grouped into sections:

```ini
[control]
alloc_model = allometric
deciduous_model = false
```

or

```ini
[files]
cfg_fname = params/NCEAS_DUKE_model_youngforest_amb.cfg
met_fname = met_data/DUKE_met_data_amb_co2.csv
out_fname = outputs/D1GDAYDUKEAMB.csv
out_param_fname = params/NCEAS_DUKE_model_simulation_amb.cfg
```

As all the model parameters are accessible via this file, these files can be quite long. Clearly it isn't necessary to list every parameter. The recommended approach is to use the [base file](example/params/base_start.cfg) and then customise whichever parameters are required via a shell script, e.g. see the python [wrapper script](example/example.py). This file just lists the parameters which needs to be changed and calls [adjust_gday_param_file.py](scripts/adjust_gday_param_file.py) to swap the parameters (listed as a python dictionary) with the default parameter. Clearly it would be trivial to write an alternative version in another language. I should highlight that I wouldn't necessarily trust the default values :).

Finally, (I need to check), but the options to vary the state and flux variables on the fly is a nice hangover from the python implementation. This functionality doesn't actually exist in the C code, instead all the state and flux variables used in the FACE intercomparisons are dumped as standard.

When I have time I will write something more extensive (ha), but information about what different variable names refer to are listed in the [header file](src/include/gday.h), which documents the different structures (i.e. control, state, params). 


## Meteorological driving file

**30-minute file:**

Variable | Description | Units
--- | --- | ---
year | |
doy  | day of year  | [0-365/6]
rain | rainfall | mm 30 min<sup>-1</sup>
par | photosynthetically active radiation | umol m<sup>-2</sup> s<sup>-1</sup>
tair | air temperature | deg C
tsoil | soil temperature | deg C
vpd | vapour pressure deficit | kPa
co2 | CO<sub>2</sub> concentration | ppm
ndep | nitrogen deposition | t ha<sup>-1</sup> 30 min<sup>-1</sup>
wind | wind speed | m<sup>-2</sup> s<sup>-1</sup>
press | atmospheric pressure | kPa


**Day file:**

Variable | Description | Units
--- | --- | ---
year | |
doy  | day of year  | [0-365/6]
tair | (daylight) air temperature | deg C
rain | rainfall | mm 30 min<sup>-1</sup>
tsoil | soil temperature | deg C
tam | morning air temperature | deg C
tpm | afternoon air temperature | deg C
tmin | minimum (day) air temperature | deg C
tmax | minimum (day) air temperature | deg C
tday | day average air temperature (24 hrs) | deg C
vpd_am | morning vapour pressure deficit | kPa
vpd_pm | afternoon vapour pressure deficit | kPa
vpd | daylight average vapour pressure deficit | kPa
co2 | CO<sub>2</sub> concentration | ppm
ndep | nitrogen deposition | t ha<sup>-1</sup> 30 min<sup>-1</sup>
wind | wind speed | m<sup>-2</sup> s<sup>-1</sup>
atmos_press | atmospheric pressure | kPa
wind_am | morning wind speed | m<sup>-2</sup> s<sup>-1</sup>
wind_pm | afternoon wind speed | m<sup>-2</sup> s<sup>-1</sup>
par | daylight photosynthetically active radiation | umol m<sup>-2</sup> s<sup>-1</sup>
par_am | morning photosynthetically active radiation | umol m<sup>-2</sup> s<sup>-1</sup>
par_am | afternoon photosynthetically active radiation | umol m<sup>-2</sup> s<sup>-1</sup>


## Example run
The [example](example) directory has two python scripts which provide an example of how one might set about running the model. [example.py](example.py) simulates the DUKE FACE experiment and [run_experiment.py](run_experiment.py) is just a wrapper script around this which produces a plot at the end comparing the data to the observations.

This example tends to break from time to time when I change various options, so please let me know if it isn't working!


## Key References
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. *Ecological Applications*, 3, 666-681.
2. Medlyn, B. E., McMurtrie, R. E., Dewar, R. C. and Jeffreys, M. P. (2000), Soil processes dominate the long-term response of forest net primary productivity to increased temperature and atmospheric CO2 concentration, *Canadian Journal of Forest Research*, 30, 873–888.
3. Sands PJ (1995) Modelling canopy production. II. From single-leaf photosynthetic parameters to daily canopy photosynthesis. Australian Journal of *Plant Physiology*, 22, 603-614.
4. Sands PJ (1996) Modelling canopy production. III. Canopy light-utilisation efficiency and its sensitivity to physiological environmental variables. *Australian Journal of Plant Physiology*, 23, 103-114.
5. de Pury, D.G.G., Farquhar, G.D. (1997) Simple scaling
of photosynthesis from leaves to canopies without the errors of big-leaf models. *Plant, Cell and Environment*, 20, 537-557.
6. Wang, Y-P. and Leuning R. (1998) A two-leaf model for canopy conductance, photosynthesis and portioning of available energy I: Model description and comparison with a multi-layered model. *Agricultural and Forest Meteorology*, 91, 89–111.

## Contacts
* [Martin De Kauwe](http://mdekauwe.github.io/).
* [Belinda Medlyn](<http://bio.mq.edu.au/people/person.php?user=bmedlyn).
