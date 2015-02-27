=============================================
GDAY model
=============================================

GDAY (Generic Decomposition And Yield) is a simple, daily time step ecosystem model that represents carbon, nitrogen, and water dynamics at the stand scale. The model is coded entirely in C without any dependancies, unless you wish to use the optimal allocation routines (see below). The wrapper files
for the example scripts and the script used to change GDAY parameter options,
are written in python. The old python version is still `online <https://github.com/mdekauwe/pygday>`_.

Key References
==============
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. *Ecological Applications*, 3, 666-681.
2. Medlyn, B. E., McMurtrie, R. E., Dewar, R. C. and Jeffreys, M. P. (2000), Soil processes dominate the long-term response of forest net primary productivity to increased temperature and atmospheric CO2 concentration, *Canadian Journal of Forest Research*, 30, 873–888.

**Note** there are many subtle changes from those original papers included in the code.

Dependancies
=============
If using the optimal allocation routines then you will need to build GDAY against the `GNU scientific libraries <http://www.gnu.org/software/gsl/>`_. Have a look at the Makefile, i.e. uncomment lines. Note you will have to set the appropriate path to the libraries on your machine.


.. contents:: :local:

Installation
=============
There is a Makefile in the src directory... ::

    $ make clean ; make 


I haven't done anything about an installation directory so you will need to move
the executable yourself.

Running the model
=================
A simple model usage can be displayed by calling GDAY as follows ::

    $ gday -u

    
Contacts
========
* `Martin De Kauwe <http://mdekauwe.github.io/>`_  (mdekauwe at gmail.com)
* `Belinda Medlyn <http://bio.mq.edu.au/people/person.php?user=bmedlyn>`_ (bmedlyn at bio.mq.edu.au).
