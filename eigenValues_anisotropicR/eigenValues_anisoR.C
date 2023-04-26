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

#include "eigenValues_anisoR.H"
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
    defineTypeNameAndDebug(eigenValues_anisoR, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        eigenValues_anisoR,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::eigenValues_anisoR::calc()
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

		volScalarField kcopy = 0.5 * tr(ReySt);

volSymmTensorField identityMat(IOobject("identityMat", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedSymmTensor("identityMat", dimensionSet(0, 0, 0, 0, 0, 0, 0), Foam::symmTensor::I));

                tempR = tempR  - (2 * identityMat * kcopy/ 3);


                forAll(Ucopy.internalField(), cellI)
                {
                        if(kcopy.internalField()[cellI] != 0)
                                tempR.ref()[cellI] = tempR.ref()[cellI] / (2 * kcopy.internalField()[cellI]);
                        else
                                tempR.ref()[cellI] = 0.0;
                }

                forAll(Ucopy.boundaryField(), patchI)
                {
                        forAll(Ucopy.boundaryField()[patchI], faceI)
                        {
                                if(kcopy.boundaryField()[patchI][faceI] != 0)
                                        tempR.boundaryFieldRef()[patchI][faceI] = (tempR.boundaryFieldRef()[patchI][faceI] / (2 * kcopy.boundaryField()[patchI][faceI]));
                                else
                                        tempR.boundaryFieldRef()[patchI][faceI] = 0.0;
                        }
                }


                volVectorField eigvalues = eigenValues(tempR);

                forAll(eigvalues.ref(), cellI)
                {
                        for(int i =0; i < 3; i++)
                        {
                                for(int j = (i+1); j < 3; j++)
                                {
                                        if(eigvalues.ref()[cellI].component(i) < eigvalues.ref()[cellI].component(j))
                                        {
                                                scalar temp = eigvalues.ref()[cellI].component(i);
                                                eigvalues.ref()[cellI].component(i) = eigvalues.ref()[cellI].component(j);
                                                eigvalues.ref()[cellI].component(j) = temp;
                                        }
                                }
                        }
                }


	        return store(resultName_, eigenValues * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::eigenValues_anisoR::eigenValues_anisoR
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

Foam::functionObjects::eigenValues_anisoR::~eigenValues_anisoR()
{}


// ************************************************************************* //
