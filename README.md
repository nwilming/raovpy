raovpy
======

Compute a repeated measures anova in python by using R.

This package allows you to compute within-subject repeated measures
analysis of variance. It uses R (via rpy2) to do so. 

The package is tightly linked to occupy's datamat. The lm_anova 
function (which also computes Greenhouse-Geisser and Huynh-Feldt
corrections) uses datamats as input. However, the aov function 
can be used without. The difference is:

    * lm_anova -> Constructs a linear model in R and then uses Anova
                from car package to compute the appropriate statistics

    * aov -> Uses R's aov function to compute the anova


The functions are somewhat checked against SPSS. However, I only
compared one dataset which gives pretty much the same result, but
that is definetly not enough.

Documentation
-------------

Only existent in the source file. Contact me (Niklas) if you have questions.


Copyright & License
-------------------

Copyright (C) 2012, Niklas Wilming
 
Licensed under GPLv2 or later, see http://www.gnu.org/licenses/gpl-2.0.txt
for license text.

