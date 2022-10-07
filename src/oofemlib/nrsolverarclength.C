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

#include "nrsolverarclength.h"
#include "verbose.h"
#include "timestep.h"
#include "mathfem.h"
// includes for ddc - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "element.h"
#include "generalboundarycondition.h"
#include "dof.h"
#include "function.h"
#include "linesearch.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"

#ifdef __PETSC_MODULE
#include "petscsolver.h"
#include "petscsparsemtrx.h"
#endif

#include <cstdio>

//#include "arclengthstd.h"

namespace oofem {


REGISTER_SparseNonLinearSystemNM( NRSolverArcLength )

    NRSolverArcLength ::NRSolverArcLength( Domain *d, EngngModel *m ) :
    NRSolver( d, m ),
    IsArcLengthBC( false )
{
    // Find out if arc length method is applied and obtain the BC pointer
    for ( size_t i = 0; i < domain->giveBcs().size(); i++ ) {
        auto &bc = domain->giveBcs()[i];
        auto salm = dynamic_cast<StandardArcLengthMethod *>( bc.get() );
        if ( salm ) {
            this->IsArcLengthBC = true;
            this->ALMtype = "standard";
            break;
        }
    }
}


NRSolverArcLength ::~NRSolverArcLength()
{
}


void NRSolverArcLength ::initializeFrom( InputRecord &ir )
{
    NRSolver ::initializeFrom( ir );
}


NM_Status
NRSolverArcLength ::solve( SparseMtrx &k, FloatArray &R, FloatArray *R0,
    FloatArray &X, FloatArray &dX, FloatArray &F,
    const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
    int &nite, TimeStep *tStep )
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X.giveSize();
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO( "NRSolver: Iteration" );
        if ( rtolf.at( 1 ) > 0.0 ) {
            OOFEM_LOG_INFO( " ForceError" );
        }
        if ( rtold.at( 1 ) > 0.0 ) {
            OOFEM_LOG_INFO( " DisplError" );
        }
        OOFEM_LOG_INFO( "\n----------------------------------------------------------------------------\n" );
    }

    l = 1.0;

    NM_Status status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = R;
    if ( R0 ) {
        RT.add( *R0 );
    }

    RRT = parallel_context->localNorm( RT );
    RRT *= RRT;

    ddX.resize( neq );
    ddX.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES


    /////////////////////////////////////
    // Find out if arc length method is applied and obtain the BC pointer
    //bool IsArcLengthBC = false;
    //std::string ALMtype;
    //ArcLengthMethod *alm;
    //for ( size_t i = 0; i < domain->giveBcs().size(); i++ ) {
    //    auto &bc = domain->giveBcs()[i];
    //    alm = dynamic_cast<ArcLengthMethod *>( bc.get() );
    //    if ( alm ) {
    //        IsArcLengthBC = true;
    //        auto salm = dynamic_cast<StandardArcLengthMethod *>( alm );
    //        if ( salm ) { // If standard ALM
    //            ALMtype = "standard";
    //        } else {
    //            // Throw exception
    //        }

    //        break;
    //    }
    //}
    StandardArcLengthMethod *salm;
    for ( size_t i = 0; i < domain->giveBcs().size(); i++ ) {
        auto &bc = domain->giveBcs()[i];
        salm     = dynamic_cast<StandardArcLengthMethod *>( bc.get() );
        if ( salm ) {
            break;
        }
    }

    // Get vector of external forces, displacement increments and lamda increments
    FloatArray Fext_in, du_in, RTmod;
    double dLam_in, LamN;

    if ( IsArcLengthBC ) {
        Fext_in.resize( neq - 1 );
        du_in.resize( neq - 1 );

        for ( int i = 1; i < neq; i++ ) {
            Fext_in.at( i ) = RT.at( i );
            du_in.at( i )   = dX.at( i );
        }
        
        dLam_in = 1.118033;
        //dLam_in = alm->give_dL();
        LamN    = X.at( neq );
        this->setArcLengthPar( Fext_in, du_in, dLam_in, true, salm);
    }
    /////////////////////////////////////

    engngModel->updateComponent( tStep, NonLinearLhs, domain );
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness( k );
    }

    nite = 0;
    for ( nite = 0;; ++nite ) {
        /////////////////////////////////////
        // extract dU and dLambda and update arc length parameters
        if ( IsArcLengthBC ) {
            if ( nite > 0 ) {
                for ( int i = 1; i < neq; i++ ) {
                    du_in.at( i ) = dX.at( i );
                }
                dLam_in = dX.at( neq );

                this->setArcLengthPar( Fext_in, du_in, dLam_in, false, salm);
            }
        }
        /////////////////////////////////////


        // Compute the residual
        engngModel->updateComponent( tStep, InternalRhs, domain );
        //rhs.beDifferenceOf( RT, F );
        ///////////////////////////
        RTmod = RT;
        // If arc length method is applied, external forces are multiplied by dLambda
        if ( IsArcLengthBC ) {
            RTmod.times( dLam_in + LamN );
            RTmod.at( neq ) = 0.0;
        }
        rhs.beDifferenceOf( RTmod, F );
        ///////////////////////////

        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement( nite, k, rhs, rlm, tStep );
        }

        // convergence check
        converged = this->checkConvergence( RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag );

        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            OOFEM_WARNING( "Divergence reached after %d iterations", nite );
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            status |= NM_Success;
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG( "Maximum number of iterations reached\n" );
            break;
        }

        if ( nite > 0 || !mCalcStiffBeforeRes ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent( tStep, NonLinearLhs, domain );
                applyConstraintsToStiffness( k );
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R.zero();
            ddX = rhs;
        } else {
            //            if ( engngModel->giveProblemScale() == macroScale ) {
            //              k.writeToFile("k.txt");
            //            }
            //k.printYourself();
            //rhs.printYourself();
            //dX.printYourself();
            linSolver->solve( k, rhs, ddX );
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
            // line search
            LineSearchNM ::LS_status LSstatus;
            double eta;
            this->giveLineSearchSolver()->solve( X, ddX, F, R, R0, prescribedEqs, 1.0, eta, LSstatus, tStep );
        } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) {
            ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
            if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                OOFEM_LOG_INFO( "Constraining increment to be %e times full increment...\n", this->constrainedNRalpha );
                ddX.times( this->constrainedNRalpha );
            }
            //this->giveConstrainedNRSolver()->solve(X, & ddX, this->forceErrVec, this->forceErrVecOld, status, tStep);
        }


        /////////////////////////////////////////

        double maxInc = 0.0;
        for ( double inc : ddX ) {
            if ( fabs( inc ) > maxInc ) {
                maxInc = fabs( inc );
            }
        }

        if ( maxInc > maxIncAllowed ) {
            if ( engngModel->giveProblemScale() == macroScale ) {
                printf( "Restricting increment. maxInc: %e\n", maxInc );
            }
            ddX.times( maxIncAllowed / maxInc );
        }

        /////////////////////////////////////////


        X.add( ddX );
        dX.add( ddX );

        if ( solutionDependentExternalForcesFlag ) {
            engngModel->updateComponent( tStep, ExternalRhs, domain );
            RT = R;
            if ( R0 ) {
                RT.add( *R0 );
            }
        }


        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();

        engngModel->giveExportModuleManager()->doOutput( tStep, true );
    }

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at( i ) ) = F.at( prescribedEqs.at( i ) ) - R0->at( prescribedEqs.at( i ) ) - R.at( prescribedEqs.at( i ) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at( i ) ) = F.at( prescribedEqs.at( i ) ) - R.at( prescribedEqs.at( i ) );
        }
    }

    this->lastReactions.resize( numberOfPrescribedDofs );

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO( "\n" );
        OOFEM_LOG_INFO( "NRSolver:     Quasi reaction table                                 \n" );
        OOFEM_LOG_INFO( "NRSolver:     Node            Dof             Displacement    Force\n" );
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R.at( prescribedEqs.at( i ) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at( i ) );
            }
            lastReactions.at( i ) = reaction;
            OOFEM_LOG_INFO( "NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at( 2 * i - 1 ), prescribedDofs.at( 2 * i ),
                X.at( prescribedEqs.at( i ) ), reaction );
        }
        OOFEM_LOG_INFO( "\n" );
    }
#endif

    return status;
}


void NRSolverArcLength::setArcLengthPar( FloatArray &Fext_in, FloatArray &du_in, double &dLam_in, 
                                        bool firstEval, StandardArcLengthMethod *salm)
{    
    //StandardArcLengthMethod *salm; // This needs to be modified, works only for standard
    //if ( ALMtype == "standard" ) {
    //    salm = dynamic_cast<StandardArcLengthMethod *>( alm );
    //} else {
    //    // Throw exception
    //}

    if ( firstEval ) {
        salm->set_Fext( Fext_in );
    }
    salm->set_du( du_in );
    salm->set_dLam( dLam_in );
    salm->setBaseParameters();
};

} // end namespace oofem
