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

Class
    Foam::functionObjects::eigenValues_anisoR

Description
    Calculates and outputs the second invariant of the velocity gradient tensor
    [1/s^2].

    \f[
        eigenValues_anisoR = 0.5(sqr(tr(\nabla U)) - tr(((\nabla U) \cdot (\nabla U))))
    \f]

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    eigenValues_anisoR.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_eigenValues_anisoR_H
#define functionObjects_eigenValues_anisoR_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class eigenValues_anisoR Declaration
\*---------------------------------------------------------------------------*/

class eigenValues_anisoR
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the eigenValues_anisoR field and return true if successful
        virtual bool calc();


//	volVectorField uMean(IOobject("uMean", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh_);// dimensionedVector("uMean", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0)));

public:

    //- Runtime type information
    TypeName("eigenValues_anisoR");
    
 
    static const word& modelName();

//     const volVectorField& Ucopy = lookupObject<volVectorField>("U");

//     volVectorField uMean = Ucopy;


    // Constructors

        //- Construct from Time and dictionary
        eigenValues_anisoR
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~eigenValues_anisoR();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
