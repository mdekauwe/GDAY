=============================================
GDAY model
=============================================

GDAY (Generic Decomposition And Yield) is a simple, daily time step ecosystem model that represents carbon, nitrogen, and water dynamics at the stand scale. The model is coded entirely in C without any dependancies. The wrapper files
for the example scripts and the script used to change GDAY parameter options,
are written in python. The old python version is still `online <https://github.com/mdekauwe/pygday>`.

Key References
==============
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. *Ecological Applications*, 3, 666-681.
2. Medlyn, B. E., McMurtrie, R. E., Dewar, R. C. and Jeffreys, M. P. (2000), Soil processes dominate the long-term response of forest net primary productivity to increased temperature and atmospheric CO2 concentration, *Canadian Journal of Forest Research*, 30, 873–888.

**Note** there are many subtle changes from those original papers included in the code.



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
