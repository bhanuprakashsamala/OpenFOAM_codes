/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "epsilon.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "turbulenceFields.H"
#include "IncompressibleTurbulenceModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(epsilon, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        epsilon,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::epsilon::calc()
{

	       // epsilon (mean) calculation with instantaneous velocity as input. 
	       // using equation 3 in Estimation of the dissipation rate of turbulent kinetic energy: A review by G Wang et.al (2020)
	       // https://www.sciencedirect.com/science/article/abs/pii/S0009250920306655

	       // velocity input
               const volVectorField& Ucopy = lookupObject<volVectorField>("U");

	       //To compute the fluctuating velocity component by first computing UMean and subtracting from instant U
               static volVectorField uMean(IOobject("uMean", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedVector("uMean", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0)));

	       static scalar counter = 1;

	       uMean = (((1/counter) * Ucopy) + ((1 - (1/counter)) * uMean));

	       volVectorField uFluc = Ucopy - uMean; 

	       //fluctuating velocity component gradient calculation
               const tmp<volTensorField> tgradu(fvc::grad(uFluc));

               const volTensorField& gradu = tgradu();

	       volScalarField xx = gradu.component(tensor::XX);
	       volScalarField xy = gradu.component(tensor::XY);
	       volScalarField xz = gradu.component(tensor::XZ);
	       volScalarField yx = gradu.component(tensor::YX);
	       volScalarField yy = gradu.component(tensor::YY);
	       volScalarField yz = gradu.component(tensor::YZ);
	       volScalarField zx = gradu.component(tensor::ZX);
	       volScalarField zy = gradu.component(tensor::ZY);
	       volScalarField zz = gradu.component(tensor::ZZ);
	
	       //kinematic viscosity input
	       const dimensionedScalar& nucopy = Ucopy.db().lookupObject<IOdictionary>("transportProperties").lookup("nu");

	       //epsilon calculation
	       volScalarField epsi = nucopy * ((2 * xx * xx) + (2 * yy * yy) + (2 * zz * zz) + (yx * yx) + (xy * xy) + (xz * xz) + (zx * zx) + (zy * zy)\
			      		+ (yz * yz) + (2 * yx * xy) + (2 * zx * xz) + (2 * zy * yz));

//	       const volScalarField& nucopy = lookupObject<volScalarField>("nu");
//	       const dimensionedScalar& nucopy = Ucopy.db().lookupObject<IOdictionary>("transportProperties").lookup("nu");

	       //computing mean epsilon
	       static volScalarField epsiMean(IOobject("epsiMean", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedScalar("epsiMean", dimensionSet(0,2,-3,0,0,0,0), Foam::scalar(0)));

	       epsiMean = (((1/counter) * epsi) + ((1 - (1/counter)) * epsiMean));

	       counter++;

	       return store(resultName_, epsiMean * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::epsilon::epsilon
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName(typeName, "U");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::epsilon::~epsilon()
{}


// ************************************************************************* //
