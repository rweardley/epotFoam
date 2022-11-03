# epotFoam

OpenFOAM electric potential formulation MHD solver.

This solver was obtained from [1] under CC-BY-NC-SA, and updated to work in OpenFOAM-10.

## Build instructions

```
cd $WM_PROJECT_USER_DIR
mkdir -p applications/solvers/electromagnetics
cd applications/solvers/electromagnetics
git clone git@github.com:rweardley/epotFoam.git
mv epotFoam epotFoamAdaptive
cd epootFoamAdaptive
git checkout epotFoamAdaptive
```
Then to compile the solver
```
wmake
```

[1] Tassone, A.: Magnetic induction and electric potential solvers for incompressible MHD flows. In Proceedings of CFD with OpenSource Software, 2016, Edited by Nilsson H. http://dx.doi.org/10.17196/OS_CFD#YEAR_2016)

