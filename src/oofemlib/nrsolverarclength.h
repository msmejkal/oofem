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


#ifndef nrsolverarclength_h
#define nrsolverarclength_h

#include <set>
#include <vector>

//#include "sparselinsystemnm.h"
//#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"
#include "nrsolver.h"

#include <memory>
#include <map>
#include "arclengthstd.h"

///@name Input fields for NRSolver
//@{
#define _IFT_NRSolverArcLength_Name "nrsolverarclength"
#define _IFT_NRSolverArcLength_albcnum "albcnum"
//@}

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements Newton-Raphson Method for arc length method
 */
class OOFEM_EXPORT NRSolverArcLength : public NRSolver
{
    bool IsArcLengthBC;
    //std::shared_ptr<ArcLengthMethod> alm;
    ArcLengthMethod *alm;
    FloatArray dXsave;

public:
    NRSolverArcLength( Domain *d, EngngModel *m );
    virtual ~NRSolverArcLength();

    // Overloaded methods:
    NM_Status solve( SparseMtrx &k, FloatArray &R, FloatArray *R0,
        FloatArray &X, FloatArray &dX, FloatArray &F,
        const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
        int &nite, TimeStep * ) override;

    void initializeFrom( InputRecord &ir ) override;
    
    const char *giveClassName() const override { return "NRSolverArcLength"; }
    virtual const char *giveInputRecordName() const override { return _IFT_NRSolverArcLength_Name; }

    void setArcLengthPar(const FloatArray &Fext_in, const FloatArray &du_in, double dLam_in, bool firstEval);

    void initiate_dXsave( int n ) { this->dXsave.resize( n ); };
};
} // end namespace oofem
#endif // nrsolver_h
