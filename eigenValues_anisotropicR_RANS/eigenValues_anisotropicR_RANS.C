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

#include "eigenValues_anisotropicR_RANS.H"
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
    defineTypeNameAndDebug(eigenValues_anisotropicR_RANS, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        eigenValues_anisotropicR_RANS,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::eigenValues_anisotropicR_RANS::calc()
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


                volScalarField C1 = (eigvalues.component(0) - eigvalues.component(1));
                volScalarField C2 = ( 2 * (eigvalues.component(1) - eigvalues.component(2)));
                volScalarField C3 = eigvalues.component(0);

                forAll(C3.ref(), cellI)
                {
                        C3.ref()[cellI] = ((3 * eigvalues.ref()[cellI].component(2)) + 1.0);
                }

                forAll(C3.boundaryFieldRef(), patchI)
                {
                        forAll(C3.boundaryFieldRef()[patchI], faceI)
                        {
                                C3.boundaryFieldRef()[patchI][faceI] = ((3 * eigvalues.boundaryFieldRef()[patchI][faceI].component(2)) + 1.0);
                        }
                }

		volScalarField C1norm = mag(C1) / sqrt(pow(C1,2) + pow(C2,2) + pow(C3,2));

                volScalarField C2norm = mag(C2) / sqrt(pow(C1,2) + pow(C2,2) + pow(C3,2));

                volScalarField C3norm = mag(C3) / sqrt(pow(C1,2) + pow(C2,2) + pow(C3,2));
		

		volScalarField etaY = (C3norm * 0.866);


		return store(resultName_, etaY * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::eigenValues_anisotropicR_RANS::eigenValues_anisotropicR_RANS
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

Foam::functionObjects::eigenValues_anisotropicR_RANS::~eigenValues_anisotropicR_RANS()
{}


// ************************************************************************* //
