/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"


#include "EigenMatrix.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Simple test ODE system
class chemistryModel
:
    public ODESystem
{

public:

    label nEqns() const
    {
        return 4;
    }

    void derivatives
    (
        const scalar x,
        const scalarField& y,
        const label li,
        scalarField& dydx
    ) const
    {
        dydx[0] = -y[1];
        dydx[1] = y[0] - (1.0/x)*y[1];
        dydx[2] = y[1] - (2.0/x)*y[2];
        dydx[3] = y[2] - (3.0/x)*y[3];
    }

    virtual void jacobian
    (
        const scalar x,
        const scalarField& y,
        const label li,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx[0] = 0.0;
        dfdx[1] = (1.0/sqr(x))*y[1];
        dfdx[2] = (2.0/sqr(x))*y[2];
        dfdx[3] = (3.0/sqr(x))*y[3];

        dfdy(0, 0) = 0.0;
        dfdy(0, 1) = -1.0;
        dfdy(0, 2) = 0.0;
        dfdy(0, 3) = 0.0;

        dfdy(1, 0) = 1.0;
        dfdy(1, 1) = -1.0/x;
        dfdy(1, 2) = 0.0;
        dfdy(1, 3) = 0.0;

        dfdy(2, 0) = 0.0;
        dfdy(2, 1) = 1.0;
        dfdy(2, 2) = -2.0/x;
        dfdy(2, 3) = 0.0;

        dfdy(3, 0) = 0.0;
        dfdy(3, 1) = 0.0;
        dfdy(3, 2) = 1.0;
        dfdy(3, 3) = -3.0/x;
    }
};

class cema
: public chemistryModel
{
protected:
    scalar cem(const scalarSquareMatrix& J) const
    {   
        // Compute Eigen values
        EigenMatrix<doubleScalar> EM(J, false);
        DiagonalMatrix<scalar> EValsRe(EM.EValsRe());
        DiagonalMatrix<scalar> EValsIm(EM.EValsIm());

        // Compute magnitude of eigen values
        scalarList EValsMag(EValsRe*EValsRe + EValsIm*EValsIm);

        // Get the indices of magnitude of eigenvalues in increasing order
        labelList order;
        sortedOrder(EValsMag, order);

        // Skip conservation modes for elements and temperature
        label nElements = 2;
        for (label i=0; i<nElements+1; ++i) 
        {
            EValsRe[order[i]] = -1E30;
        }

        // return cem
        return max(EValsRe);
    }

public:
    
    virtual void jacobian
    (
        const scalar x,
        const scalarField& y,
        const label li,
        scalarField& dfdx,
        scalarSquareMatrix& J
    ) const
    {
        chemistryModel::jacobian(x, y, li, dfdx, J);
        cem(J);
    }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Create the ODE system
    cema ode;

    // Initialize nEq of the system
    const label nEq = ode.nEqns();

    // Initialize dummy state vector y at some x
    scalar x(1.0);
    scalarField y(nEq, 1.0);

    // Compute the jacobian matrix
    scalarField dfdx(nEq);
    scalarSquareMatrix dfdy(nEq);
    ode.jacobian(x, y, 0, dfdx, dfdy);
    
    return 0;
}


// ************************************************************************* //
