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

#include "arclength.h"
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
REGISTER_BoundaryCondition( ArcLengthMethod );

ArcLengthMethod ::ArcLengthMethod( int n, Domain *d ) :
    ActiveBoundaryCondition( n, d ),
    md( new Node( 0, domain ) ),
    w(0.0),
    g(0.0)
{
    // this allocates a new equation related to arc length method
    this->md->appendDof( new MasterDof( this->md.get(), (DofIDItem)( d->giveNextFreeDofID() ) ) );
}


ArcLengthMethod ::~ArcLengthMethod()
{
}


void ArcLengthMethod ::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition ::initializeFrom( ir );
    //IR_GIVE_FIELD( ir, dL, _IFT_ArcLengthMethod_dL );
}


void ArcLengthMethod ::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s,
    const UnknownNumberingScheme &c_s, double scale,
    void *lock )
{
    int nEq = r_s.giveDofEquationNumber( *md->begin() ); // number of last equation
    IntArray nEqArr( 1 );
    nEqArr.at( 1 ) = nEq;

    IntArray locArr( nEq - 1 ); // location array
    for ( int i = 1; i < nEq; i++ ) {
        locArr.at( i ) = i;
    }

    FloatMatrix wMat( 1,1 );
    wMat.at( 1,1 ) = w;

    // Create FloatMatrix from FloatArray
    FloatArray FextNegArr = Fext;
    FextNegArr.times( -1 );
    FloatMatrix FextNeg( FextNegArr, 0 );

    answer.assemble( nEqArr, locArr, this->h );
    answer.assemble( locArr, nEqArr, FextNeg );
    answer.assemble( nEqArr, nEqArr, wMat );
}

void ArcLengthMethod ::assembleVector( FloatArray &answer, TimeStep *tStep,
    CharType type, ValueModeType mode,
    const UnknownNumberingScheme &s,
    FloatArray *eNorms,
    void *lock )
{
    int nEq = s.giveDofEquationNumber( *md->begin() ); // number of last equation
    IntArray nEqArr( 1 );
    nEqArr.at( 1 ) = nEq; 

    // rhs
    FloatArray gNegArr( 1 );
    gNegArr.at( 1 ) = g;

    answer.assemble( gNegArr, nEqArr );
}

} // namespace oofem
