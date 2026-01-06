## Description

OpenFOAM solver for mechanistic nucleate pool boiling, capturing interface-localized evaporation/condensation, stochastic nucleation, contact angle boundary condition, and conjugate heat transfer. Accurately predicts bubble growth, detachment, coalescence, and vapor-cluster formation, validated against FC-72/FC-77 experiments.



## Boundary condition instructions

Replace the derivedFvPatchFields directory located at
src/regionModels/thermalBaffleModels/derivedFvPatchFields with the provided version.
After replacement, run wmake to recompile.

This enables the use of the boundary condition thermalBaffleTop, which allows extension of the bottom boundary through the thermal baffle model.

This allows to solve conjugate heat transfer (both solid and fluid domain)



## Solver execution

Go to folder BoilingCondensationFOAM\_solver

Execute the ./compile\_command to create a compiled solver



## License
This project is licensed under the **GNU GPL v3.0** — see the `LICENSE` file for details.
