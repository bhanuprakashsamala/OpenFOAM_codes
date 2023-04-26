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

#include "R_usingVelocity.H"
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
    defineTypeNameAndDebug(R_usingVelocity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        R_usingVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::R_usingVelocity::calc()
{

               const volVectorField& Ucopy = lookupObject<volVectorField>("U");

               static volVectorField uMean(IOobject("uMean", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedVector("uMean", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0)));

	       static scalar counter = 1;

	      static volSymmTensorField ReySt(IOobject("ReySt", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedSymmTensor("ReySt", dimensionSet(0,2,-2,0,0,0,0), symmTensor(0,0,0,0,0,0)));

	       ReySt += sqr(uMean);

	       uMean = (((1/counter) * Ucopy) + ((1 - (1/counter)) * uMean));

	       scalar beta = 1/counter;

               ReySt = ((1 - beta ) * ReySt) + (beta * sqr(Ucopy)) - (sqr(uMean));

	       counter++;
	
	       volSymmTensorField tempR = ReySt; 

	       return store(resultName_, tempR * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::R_usingVelocity::R_usingVelocity
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

Foam::functionObjects::R_usingVelocity::~R_usingVelocity()
{}


// ************************************************************************* //
