/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    epotFoam

Description
    Transient solver for incompressible, laminar MHD flow of Newtonian fluids
    employing the electric potential formulation. The Lorentz force term is
    treated in a NON-CONSERVATIVE way, but a consistent and conservative scheme
    for the computation of the current density according to Ni et al. is used to
    ensure conservation of momentum and current density.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Lorentz force term initialization 
    volVectorField lorentz = sigma * (-fvc::grad(PotE) ^ B0) + sigma * ((U ^ B0) ^ B0);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
          - (1.0/rho) * lorentz //Lorentz term
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            // field must satisfy continuity equation
            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn // Pressure Equation coefficient matrix
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                // reference value of pressure
                pEqn.setReference(pRefCell, pRefValue);

                // solve pressure equation according to solution dictionary, check if final iteration
                // pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        //Interpolating cross product u x B over mesh faces 
        surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();

        //Poisson equation for electric potential 
        fvScalarMatrix PotEEqn
        (
            fvm::laplacian(PotE)
            ==
            fvc::div(psiub)
        );

        //Reference potential 
        PotEEqn.setReference(PotERefCell, PotERefValue);

        //Solving Poisson equation
        PotEEqn.solve();

        //Computation of current density at cell faces
        surfaceScalarField jn = -(fvc::snGrad(PotE) * mesh.magSf()) + psiub;
        
        //Current density at face center
        surfaceVectorField jnv = jn * mesh.Cf();
        
        //Interpolation of current density at cell center 
        volVectorField jfinal = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());
        
        //Update current density distribution and boundary condition 
        jfinal.correctBoundaryConditions();

        //Lorentz force computation 
        lorentz = sigma* (jfinal ^ B0);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
