## Information 

A new field variable Cc_mod is introduced to modulate condensation strength. The effective accommodation coefficient used in the condensation model is:

𝐶oeffC_eff = Cc_mod×CoeffC

To enforce near-complete vapor return, set Cc_mod to a very large value, or define a localized region near the top wall with very high Cc_mod, so that any vapor reaching the top boundary condenses back into the domain

## Execution

<<<<<<< HEAD
Execute the solver : BoilingCondensationFOAM
=======
Execute the solver : BoilingCondensationFOAM


## Changes in constant folder

Always keep 
Cc 0;
Cv 0;

//ignore these three lines, these are for cavitation only please 
n             1.6e+10; // required for cavitation only, so not here. 
dNuc          2.0e-06; // required for cavitation only, so not here. 
PSat	   1e5; 
>>>>>>> 4a4ddcf (Initial Upload)
