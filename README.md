# epotFoam

OpenFOAM electric potential formulation MHD solver.

This solver was obtained from [1] under CC-BY-NC-SA, and updated to work in OpenFOAM-10.

This version has been further modified to treat the imposed magnetic field B0 as a volVectorField set in <case>/0/B0, rather than as a constant in <case>/constant/electromagneticProperties. 
This allows for the imposed magnetic field to be non-uniform in space, by writing a codeStream.

## Build instructions

```
cd $WM_PROJECT_USER_DIR
mkdir -p applications/solvers/electromagnetics
cd applications/solvers/electromagnetics
git clone git@github.com:rweardley/epotFoam.git
mv epotFoam epotFoamNonuniform
cd epootFoamNonuniform
git checkout epotFoamNonuniform
```
Then to compile the solver
```
wmake
```

[1] Tassone, A.: Magnetic induction and electric potential solvers for incompressible MHD flows. In Proceedings of CFD with OpenSource Software, 2016, Edited by Nilsson H. http://dx.doi.org/10.17196/OS_CFD#YEAR_2016)

