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

#ifndef isolinearelasticmaterialtemperature_h
#define isolinearelasticmaterialtemperature_h

#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "scalarfunction.h" // Modification

//#include "qcmaterialextensioninterface.h"

//#include "MixedPressure/mixedpressurematerialextensioninterface.h"

///@name Input fields for IsotropicLinearElasticMaterialTemperature
//@{
#define _IFT_IsotropicLinearElasticMaterialTemperature_Name "isoletemp"
#define _IFT_IsotropicLinearElasticMaterialTemperature_e "e"
#define _IFT_IsotropicLinearElasticMaterialTemperature_n "n"
#define _IFT_IsotropicLinearElasticMaterialTemperature_talpha "talpha"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements an isotropic linear elastic material with parameters dependent on 
 * temperature in a finite element problem.
 * For large deformation analysis it becomes the St. Venant-Kirchoff hyperelasticity model.
 *
 * Tasks:
 * - Returning standard material stiffness matrix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - methods give2dPlaneStressMtrx, givePlaneStrainMtrx, give1dStressMtrx are
 *   introduced since form of this matrices is well known, and for
 *   faster response mainly in linear elastic problems.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class IsotropicLinearElasticMaterialTemperature : public LinearElasticMaterial
{
protected:
    /// Young's modulus.
    ScalarFunction E = 0;
    /// Poisson's ratio.
    ScalarFunction nu = 0;
    /// Alpha
    double a = 0;

public:
    /**
     * Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param n material model number in domain
     * @param d domain which receiver belongs to
     */
    IsotropicLinearElasticMaterialTemperature( int n, Domain *d );
    /**
     * Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param n Material model number in domain.
     * @param d Domain which receiver belongs to.
     * @param E Young modulus.
     * @param nu Poisson ratio.
     */
    IsotropicLinearElasticMaterialTemperature( int n, Domain *d, double E, double nu );

    const char *giveClassName() const override { return "IsotropicLinearElasticMaterialTemperature"; }
    const char *giveInputRecordName() const override { return _IFT_IsotropicLinearElasticMaterialTemperature_Name; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The E modulus (keyword "E"), Poisson ratio ("nu") and coefficient of thermal dilatation
     * alpha ("talpha") are read. The parent class instanciateFrom method is called.
     */
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    //double give(int aProperty, GaussPoint *gp) const override;

    /// Returns Young's modulus.
    double giveYoungsModulus( GaussPoint *gp, TimeStep *tStep ) const { return eval_E(gp, tStep ); }

    /// Returns Poisson's ratio.
    double givePoissonsRatio( GaussPoint *gp, TimeStep *tStep ) const { return eval_nu( gp, tStep ); }

    /// Returns the shear elastic modulus @f$ G = \frac{E}{2(1+\nu)} @f$.
    //double giveShearModulus( ) const override { return G; } // Makes not much sense

    /// Returns the bulk elastic modulus @f$ K = \frac{E}{3(1-2\nu)} @f$.
    double giveBulkModulus( GaussPoint *gp, TimeStep *tStep ) const { return eval_E( gp, tStep ) / ( 3. * ( 1. - 2. * eval_nu( gp, tStep ) ) ); }

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode, GaussPoint * gp,
                                               TimeStep * tStep)  const override;


    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode, GaussPoint * gp,
                                               TimeStep * tStep) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode, GaussPoint * gp,
                                            TimeStep * tStep) const override;

    FloatMatrixF<6, 6> give3dMaterialStiffnessMatrix( MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep ) const override;


    FloatArrayF<6> giveRealStressVector_3d( const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep ) const override;

    FloatArrayF<3> giveRealStressVector_PlaneStress( const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep ) const override;

    FloatArrayF<1> giveRealStressVector_1d( const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep ) const override;

    // Plane Strain case is missing!!! 


    /**
     * Computes bulk modulus from given Young's modulus and Poisson's ratio.
     * @param young Young's modulus (@f$ E @f$).
     * @param nu Poisson's ratio (@f$ \nu @f$).
     * @return Bulk modulus (@f$ K = E/(3*(1-2*nu) @f$).
     */
    //static double computeBulkModulusFromYoungAndPoisson(double young, double nu)
    //{
    //    return young / ( 3. * ( 1. - 2. * nu ) );
    //}

    /**
     * Computes shear modulus from given Young's modulus and Poisson's ratio.
     * @param young Young's modulus (@f$ E @f$).
     * @param nu Poisson's ratio (@f$ \nu @f$).
     * @return Shear modulus (@f$ G = \frac{E}{2 (1+\nu)} @f$).
     */
    /*static double computeShearModulusFromYoungAndPoisson(double young, double nu)
    {
        return young / ( 2. * ( 1. + nu ) );
    }*/



    double giveTemperature( GaussPoint *gp, TimeStep *tStep ) const;

    double eval_E( GaussPoint *gp, TimeStep *tStep ) const;

    double eval_nu( GaussPoint *gp, TimeStep *tStep ) const;

    double compute_G( GaussPoint *gp, TimeStep *tStep ) const {
        return eval_E( gp, tStep ) / ( 2.0 * ( 1. + eval_nu( gp, tStep ) ) );
    };
    
};
} // end namespace oofem
#endif // isolinearelasticmaterial_h
