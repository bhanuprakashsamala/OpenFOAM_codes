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

#include "eigenVecsSort.H"
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
    defineTypeNameAndDebug(eigenVecsSort, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        eigenVecsSort,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::eigenVecsSort::calc()
{
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");

		const volScalarField& kcopy = Ucopy.db().lookupObject<volScalarField>("k");

		typedef incompressible::turbulenceModel icoModel;

		tmp<volSymmTensorField> Rcopy;

		tmp<volSymmTensorField> Rtemp;

		if(mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
		{
			const icoModel& model = mesh_.lookupObject<icoModel>(turbulenceModel::propertiesName);
			Rcopy = model.devReff();
			Rtemp = model.R();
		}
		else
		{
			FatalErrorInFunction << "Unable to find turbulence model in the database" << exit(FatalError);
		}

		volSymmTensorField tempR = Rtemp();

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

		
		volTensorField eigvecs(IOobject("eigvecs", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedTensor("eigvecs", dimensionSet(0, 0, 0, 0, 0, 0, 0), Foam::tensor::I)); 

		forAll(eigvecs.ref(), cellI)
		{
			eigvecs.ref()[cellI] = eigenVectors(tempR.ref()[cellI]);
		}

		forAll(eigvecs.boundaryFieldRef(), patchI)
		{
			forAll(eigvecs.boundaryFieldRef()[patchI], faceI)
			{
				eigvecs.boundaryFieldRef()[patchI][faceI] = eigenVectors(tempR.boundaryFieldRef()[patchI][faceI]);
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
						for(int k = 0; k < 3; k++)
						{
							scalar temp2 = eigvecs.ref()[cellI].component((3* i) + k);
							eigvecs.ref()[cellI].component((3*i) + k) = eigvecs.ref()[cellI].component((3*j) + k);
							eigvecs.ref()[cellI].component((3*j) + k) = temp2;
						}
                                        }
                                }
                        }
                }
		return store(resultName_, eigvecs * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::eigenVecsSort::eigenVecsSort
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

Foam::functionObjects::eigenVecsSort::~eigenVecsSort()
{}


// ************************************************************************* //
