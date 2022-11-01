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

#ifndef arclengthstd_h
#define arclengthstd_h

#include "arclength.h"
#include "floatarray.h"
#include "intarray.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"
#include "floatmatrix.h"

#include <memory>

#define _IFT_StandardArcLengthMethod_Name "arclenghtstd"

///@name Input fields for active boundary condition
//@{
#define _IFT_StandardArcLengthMethod_dL "dl"
#define _IFT_StandardArcLengthMethod_Psi "psi"
#define _IFT_StandardArcLengthMethod_dLam0 "dlam0"
//@}

namespace oofem {
/**
*Class implementing the standard arclenght method
 */
class OOFEM_EXPORT StandardArcLengthMethod : public ArcLengthMethod
{
protected:
    double dL;
    double Psi;
    FloatArray du;
    double dLam;
    double dLam0;

public:
    StandardArcLengthMethod( int n, Domain *d );
    /// Destructor.
    virtual ~StandardArcLengthMethod();

    void initializeFrom( InputRecord &ir ) override;
    const char *giveInputRecordName() const override { return _IFT_StandardArcLengthMethod_Name; }

    const char *giveClassName() const override { return "StandardArcLengthMethod"; };


    void compute_g() override;
    void compute_H() override;
    void compute_w() override;

    void set_du( const FloatArray &du_in ) override { this->du = du_in; };
    void set_dLam( const double &dLam_in ) override { this->dLam = dLam_in; };


    double give_dL() { return dL; };
    double givePsi() { return Psi; };
    double give_dLam0() override { return dLam0; };

    double computIntialGuess( FloatArray &dX, const FloatArray &d, const FloatArray &dXsave ) override;

};
} // namespace oofem
#endif // LinearConstraintBC_h
