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

\*---------------------------------------------------------------------------*/

#include "SVD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Initialize Rectangular Matrix
    RectangularMatrix<scalar> A(2, 2);
    A(0, 0) = 3.0;
    A(0, 1) = 2.0;
    A(1, 0) = 2.0;
    A(1, 1) = 3.0;

    // Print the Matrix
    Info << A << endl;

    // Compute SVD
    SVD svd(A);

    // Print singular values
    Info << svd.U() << endl;
    Info << svd.V() << endl;
    Info << svd.S() << endl;

    Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
