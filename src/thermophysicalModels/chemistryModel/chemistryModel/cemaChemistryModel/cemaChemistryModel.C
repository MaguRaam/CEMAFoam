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

\*---------------------------------------------------------------------------*/

#include "cemaChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::cemaChemistryModel<ThermoType>::cemaChemistryModel
(
    const fluidMulticomponentThermo& thermo
)
:
    chemistryModel<ThermoType>(thermo),
    nElements_(this->template lookup<label>("nElements"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::cemaChemistryModel<ThermoType>::~cemaChemistryModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::cemaChemistryModel<ThermoType>::cema
(
    const scalarSquareMatrix& J
) const
{
    // Compute Eigen values
    const EigenMatrix<scalar> EM(J, false);
    DiagonalMatrix<scalar> EValsRe(EM.EValsRe());
    const DiagonalMatrix<scalar> EValsIm(EM.EValsIm());

    // Compute magnitude square of eigen values
    const scalarList EValsMag(EValsRe*EValsRe + EValsIm*EValsIm);

    // Get the indices of magnitude of eigenvalues in increasing order
    labelList order;
    sortedOrder(EValsMag, order);

    // Skip conservation modes for elements and temperature
    const scalar smallestEVal = -vGreat;
    
    for (label i = 0; i < nElements_ + 1; ++i)
    {
        EValsRe[order[i]] = smallestEVal;
    }

    // return cem
    return max(EValsRe);
}


template<class ThermoType>
void Foam::cemaChemistryModel<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt,
    scalarSquareMatrix& J
) const
{
    chemistryModel<ThermoType>::jacobian(t, YTp, li, dYTpdt, J);
    Info << "CEMA = " << cema(J) << endl;
    Info << nElements_ << endl;
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
