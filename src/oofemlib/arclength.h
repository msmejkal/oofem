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

#ifndef arclength_h
#define arclength_h

#include "activebc.h"
#include "floatarray.h"
#include "intarray.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"
#include "floatmatrix.h"

#include <memory>

#define _IFT_ArcLengthMethod_Name "arclenght"

///@name Input fields for active boundary condition
//@{



//@}

namespace oofem {
/**
*Class implementing the arclenght method
 */
class OOFEM_EXPORT ArcLengthMethod : public ActiveBoundaryCondition
{
protected:
    std ::unique_ptr<DofManager> md;


    FloatMatrix h;
    double w;
    double g;
    FloatArray Fext;

public:
    ArcLengthMethod( int n, Domain *d );
    /// Destructor.
    virtual ~ArcLengthMethod();

    void initializeFrom( InputRecord &ir ) override;
    const char *giveInputRecordName() const override { return _IFT_ArcLengthMethod_Name; }
    void assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type, const UnknownNumberingScheme &r_s,
        const UnknownNumberingScheme &c_s, double scale = 1.0,
        void *lock = nullptr ) override;
    void assembleVector( FloatArray &answer, TimeStep *tStep,
        CharType type, ValueModeType mode,
        const UnknownNumberingScheme &s, FloatArray *eNorms = nullptr,
        void *lock = nullptr ) override;


    void saveContext( DataStream &stream, ContextMode mode ) override {};
    void restoreContext( DataStream &stream, ContextMode mode ) override {};

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager( int i ) override { return this->md.get(); }


    const char *giveClassName() const override { return "ArcLengthMethod"; };

    virtual void compute_g(){};
    virtual void compute_H(){};
    virtual void compute_w(){};  

    virtual void set_du( FloatArray &du_in ) = 0;
    virtual void set_dLam( double &dLam_in ) = 0;

    void set_Fext( FloatArray &Fext_in ) { this->Fext = Fext_in; };


};
} // namespace oofem
#endif // LinearConstraintBC_h
