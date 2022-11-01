/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "arclengthstd.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "node.h"
#include "domain.h"


namespace oofem {
REGISTER_BoundaryCondition( StandardArcLengthMethod );

StandardArcLengthMethod ::StandardArcLengthMethod( int n, Domain *d ) :
    ArcLengthMethod( n, d )
{
}

StandardArcLengthMethod ::~StandardArcLengthMethod()
{
}


void StandardArcLengthMethod ::initializeFrom( InputRecord &ir )
{
    ArcLengthMethod ::initializeFrom( ir );
    IR_GIVE_FIELD( ir, dL, _IFT_StandardArcLengthMethod_dL );
    IR_GIVE_FIELD( ir, Psi, _IFT_StandardArcLengthMethod_Psi );
    IR_GIVE_FIELD( ir, dLam0, _IFT_StandardArcLengthMethod_dLam0);
}

void StandardArcLengthMethod ::compute_g()
{
    this->g = this->du.computeSquaredNorm() + this->dLam * this->dLam * this->Psi * this->Psi * this->Fext.computeSquaredNorm() - this->dL * this->dL;
}

void StandardArcLengthMethod ::compute_H()
{
    this->h = FloatMatrix( 2 * this->du, 1 );
}

void StandardArcLengthMethod ::compute_w()
{
    this->w = 2 * this->dLam * this->Psi * this->Psi * this->Fext.computeSquaredNorm();
}

double StandardArcLengthMethod::computIntialGuess(FloatArray& dX, const FloatArray& d, const FloatArray& dXsave) {
    int neq = dX.giveSize();

    double dLam_in1 = this->dL / sqrt( d.dotProduct( d ) + this->Psi * this->Psi * this->Fext.dotProduct( this->Fext ) );
    double dLam_in2 = -dLam_in1;

    FloatArray duProv1 = d;
    FloatArray duProv2 = d;
    duProv1.times( dLam_in1 );
    duProv2.times( dLam_in2 );
    duProv1.at( neq ) = dLam_in1;
    duProv2.at( neq ) = dLam_in2;

    double DP1 = dXsave.dotProduct( duProv1 );
    double DP2 = dXsave.dotProduct( duProv2 );

    double dLam_in;
    if ( DP1 >= DP2 ) {
        dLam_in = dLam_in1;
        dX      = duProv1;
    } else if ( DP1 < DP2 ) {
        dLam_in = dLam_in2;
        dX      = duProv2;
    }

    return dLam_in;
}

} // namespace oofem
