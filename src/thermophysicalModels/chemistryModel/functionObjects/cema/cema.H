/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::cema

Description
    This function object...

    Example of function object specification:
    \verbatim
    cema1
    {
        type           cema;
        libs ("libcemaFunctionObject.so");
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: cema | yes     |
    \endtable

SourceFiles
    cema.C

\*---------------------------------------------------------------------------*/

#ifndef cema_H
#define cema_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cema Declaration
\*---------------------------------------------------------------------------*/

class cema
:
    public fvMeshFunctionObject
{
    // Private Data

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("cema");


    // Constructors

        //- Construct from Time and dictionary
        cema
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cema(const cema&) = delete;


    //- Destructor
    virtual ~cema();


    // Member Functions

        //- Read the cema data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the cema
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cema&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
