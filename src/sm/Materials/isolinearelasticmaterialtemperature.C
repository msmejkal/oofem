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

#include "linearelasticmaterial.h"
#include "isolinearelasticmaterialtemperature.h"
#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"
#include "fieldmanager.h" // Modification
#include "sm/Elements/structuralelement.h" // Modification
#include "engngm.h" // Modification
#include "function.h" // Modification

namespace oofem {
REGISTER_Material( IsotropicLinearElasticMaterialTemperature );

IsotropicLinearElasticMaterialTemperature ::IsotropicLinearElasticMaterialTemperature( int n, Domain *d ) :
    LinearElasticMaterial( n, d )
{
}

IsotropicLinearElasticMaterialTemperature ::IsotropicLinearElasticMaterialTemperature( int n, Domain *d,
    double _E, double _nu ) :
    LinearElasticMaterial( n, d ),
    E( _E ),
    nu( _nu ),
    a( 0. )
{
    this->alpha = {
        a, a, a, 0., 0., 0.
    };
}

void IsotropicLinearElasticMaterialTemperature ::initializeFrom( InputRecord &ir)
{
    LinearElasticMaterial ::initializeFrom( ir );

    IR_GIVE_FIELD( ir, E, _IFT_IsotropicLinearElasticMaterialTemperature_e );
    IR_GIVE_FIELD( ir, nu, _IFT_IsotropicLinearElasticMaterialTemperature_n );
    IR_GIVE_FIELD( ir, a, _IFT_IsotropicLinearElasticMaterialTemperature_talpha );

    this->alpha = {
        a, a, a, 0., 0., 0.
    };
}


void IsotropicLinearElasticMaterialTemperature ::giveInputRecord( DynamicInputRecord &input )
{
    this->LinearElasticMaterial ::giveInputRecord( input );
    StructuralMaterial ::giveInputRecord( input );

    input.setField( this->E, _IFT_IsotropicLinearElasticMaterialTemperature_e );
    input.setField( this->nu, _IFT_IsotropicLinearElasticMaterialTemperature_n );
    input.setField( this->a, _IFT_IsotropicLinearElasticMaterialTemperature_talpha );
}


void IsotropicLinearElasticMaterialTemperature ::saveContext( DataStream &stream, ContextMode mode )
{
    LinearElasticMaterial ::saveContext( stream, mode );

    //if ( ( mode & CM_Definition ) ) {
    //    if ( !stream.write( E ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //    if ( !stream.write( nu ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //    if ( !stream.write( a ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //}
}

void IsotropicLinearElasticMaterialTemperature ::restoreContext( DataStream &stream, ContextMode mode )
{
    LinearElasticMaterial ::restoreContext( stream, mode );

    //if ( mode & CM_Definition ) {
    //    if ( !stream.read( E ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //    if ( !stream.read( nu ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //    if ( !stream.read( a ) ) {
    //        THROW_CIOERR( CIO_IOERR );
    //    }
    //}
}

//double
//IsotropicLinearElasticMaterialTemperature ::give( int aProperty, GaussPoint *gp ) const
//{
//    if ( ( aProperty == NYxy ) || ( aProperty == NYxz ) || ( aProperty == NYyz ) ) {
//        return nu;
//    }
//
//    if ( ( aProperty == 'G' ) || ( aProperty == Gyz ) || ( aProperty == Gxz ) || ( aProperty == Gxy ) ) {
//        return G;
//    }
//
//    if ( ( aProperty == 'E' ) || ( aProperty == Ex ) || ( aProperty == Ey ) || ( aProperty == Ez ) ) {
//        return E;
//    }
//
//    if ( ( aProperty == 'n' ) || ( aProperty == NYzx ) || ( aProperty == NYzy ) || ( aProperty == NYyx ) ) {
//        return nu;
//    }
//
//    return this->Material ::give( aProperty, gp );
//}


FloatMatrixF<3, 3>
IsotropicLinearElasticMaterialTemperature ::givePlaneStressStiffMtrx( MatResponseMode mode,
                                                                      GaussPoint *gp,
                                                                      TimeStep *tStep ) const
{
    double K      = eval_E( gp, tStep ) / ( 3.0 * ( 1. - 2. * eval_nu( gp, tStep ) ) );
    FloatMatrixF<6, 6> tangentLoc = 2 * compute_G( gp, tStep ) * I_dev6 + K * I6_I6; // 3D stiffness matrix

    auto c = inv( tangentLoc );
    FloatMatrixF<3, 3> reduced = {
        c( 0, 0 ), c( 0, 1 ), c( 0, 5 ),
        c( 1, 0 ), c( 1, 1 ), c( 1, 5 ),
        c( 5, 0 ), c( 5, 1 ), c( 5, 5 ),
    };
    FloatMatrixF<3, 3>  tangentPlaneStressLoc = inv( reduced );


    if ( ( tStep->giveIntrinsicTime() < this->castingTime ) ) {
        return tangentPlaneStressLoc * ( 1. - this->preCastStiffnessReduction );
    } else {
        return tangentPlaneStressLoc;
    }
}


FloatMatrixF<4, 4>
IsotropicLinearElasticMaterialTemperature ::givePlaneStrainStiffMtrx( MatResponseMode mode,
                                                                      GaussPoint *gp,
                                                                      TimeStep *tStep ) const
{
    double K = eval_E( gp, tStep ) / ( 3.0 * ( 1. - 2. * eval_nu( gp, tStep ) ) );
    FloatMatrixF<6, 6> tangentLoc = 2 * compute_G( gp, tStep ) * I_dev6 + K * I6_I6; // 3D stiffness matrix
    FloatMatrixF<4, 4> tangentPlaneStrainLoc = {
        tangentLoc( 0, 0 ), tangentLoc( 1, 0 ), tangentLoc( 2, 0 ), tangentLoc( 3, 0 ),
        tangentLoc( 0, 1 ), tangentLoc( 1, 1 ), tangentLoc( 2, 1 ), tangentLoc( 3, 1 ),
        tangentLoc( 0, 2 ), tangentLoc( 1, 2 ), tangentLoc( 2, 2 ), tangentLoc( 3, 2 ),
        tangentLoc( 0, 3 ), tangentLoc( 1, 3 ), tangentLoc( 2, 3 ), tangentLoc( 3, 3 ),
    };

    if ( ( tStep->giveIntrinsicTime() < this->castingTime ) ) {
        return tangentPlaneStrainLoc * ( 1. - this->preCastStiffnessReduction );
    } else {
        return tangentPlaneStrainLoc;
    }
}


FloatMatrixF<1, 1>
IsotropicLinearElasticMaterialTemperature ::give1dStressStiffMtrx( MatResponseMode mode,
                                                                   GaussPoint *gp,
                                                                   TimeStep *tStep ) const
{
    double e = this->eval_E(gp, tStep);
    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        e *= 1. - this->preCastStiffnessReduction;
    }
    return { e };
}


FloatMatrixF<6, 6>
IsotropicLinearElasticMaterialTemperature ::give3dMaterialStiffnessMatrix( MatResponseMode mode,
                                                                           GaussPoint *gp,
                                                                           TimeStep *tStep ) const
{
    double K = eval_E( gp, tStep ) / ( 3.0 * ( 1. - 2. * eval_nu( gp, tStep ) ) );
    FloatMatrixF<6, 6> tangentLoc = 2 * compute_G( gp, tStep ) * I_dev6 + K * I6_I6; // 3D stiffness matrix

    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        return tangentLoc * ( 1. - this->preCastStiffnessReduction );
    } else {
        return tangentLoc;
    }
}

FloatArrayF<6>
IsotropicLinearElasticMaterialTemperature ::giveRealStressVector_3d( const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    auto d = this->give3dMaterialStiffnessMatrix( TangentStiffness, gp, tStep );

    FloatArrayF<6> stress;

    auto thermalStrain   = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Incremental );

    auto strainIncrement = strain - thermalStrain - FloatArrayF<6>( ( status->giveStrainVector() ) );

    stress = dot( d, strainIncrement ) + status->giveStressVector();

    // update gp
    status->letTempStrainVectorBe( strain );
    status->letTempStressVectorBe( stress );
    return stress;
}

FloatArrayF<3>
IsotropicLinearElasticMaterialTemperature ::giveRealStressVector_PlaneStress( const FloatArrayF<3> &reducedStrain, 
                                                                                    GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    auto d = this->givePlaneStressStiffMtrx( TangentStiffness, gp, tStep );

    FloatArray stress, stressTest, stressDiff, stressOld;
    FloatArray strainVector, strainVectorTest, strainOldTemp;

//    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
//        this->giveStressDependentPartOfStrainVector( strainVector, gp, reducedStrain, tStep, VM_Total );
//        stress.beProductOf( d, strainVector );
//
/////////////////////////////////////////////////////////////////////
//        // Test if stresses for both formulations differ when the parameters are independent of temperature
//        FloatArray gcoords;
//
//        (gp->giveElement() )->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates()) ;
//
//        this->giveStressDependentPartOfStrainVector( strainVectorTest, gp, reducedStrain, tStep, VM_Incremental );
//
//        auto strainOld = status->giveStrainVector();
//
//        auto strainIncrementTest = strainVectorTest - strainOld;
//
//        auto stressOld = status->giveStressVector();
//
//        stressTest.beProductOf( d, strainIncrementTest );
//
//        stressTest.add( stressOld );
//
//        stressDiff = stressTest - stress;
//
//        double DiffNormSq = stressDiff.dotProduct( stressDiff );
//        if ( DiffNormSq > 1e-6) {
//            double x = 0;
//        }
/////////////////////////////////////////////////////////////////////
//
//    } else { // changes in material stiffness ->> incremental formulation
//        this->giveStressDependentPartOfStrainVector( strainVector, gp, reducedStrain, tStep, VM_Incremental );
//        auto strainIncrement = strainVector - status->giveStrainVector();
//
//        stress.beProductOf( d, strainIncrement );
//        stress.add( status->giveStressVector() );
//    }
// 
    // Incremental formulation
    this->giveStressDependentPartOfStrainVector( strainVector, gp, reducedStrain, tStep, VM_Incremental );
    FloatArray redStrain;
    StructuralMaterial::giveReducedSymVectorForm( redStrain, status->giveStrainVector(), _PlaneStress );
    auto strainIncrement = strainVector - redStrain;

    stress.beProductOf( d, strainIncrement );
    FloatArray redStress;
    StructuralMaterial::giveReducedSymVectorForm( redStress, status->giveStressVector(), _PlaneStress );
    stress.add( redStress );


    // update gp
    status->letTempStrainVectorBe( reducedStrain );
    status->letTempStressVectorBe( stress );
    return stress;
}

FloatArrayF<1>
IsotropicLinearElasticMaterialTemperature ::giveRealStressVector_1d( const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    auto d      = this->give1dStressStiffMtrx( TangentStiffness, gp, tStep );

    FloatArray answer;
    FloatArray strainVector;

    // changes in material stiffness ->> incremental formulation
    this->giveStressDependentPartOfStrainVector( strainVector, gp, reducedStrain, tStep, VM_Incremental );
    auto strainIncrement = strainVector - status->giveStrainVector();

    answer.beProductOf( d, strainIncrement );
    answer.add( status->giveStressVector() );


    // update gp
    status->letTempStrainVectorBe( reducedStrain );
    status->letTempStressVectorBe( answer );
    return answer;
}

double IsotropicLinearElasticMaterialTemperature::giveTemperature( GaussPoint *gp, TimeStep *tStep ) const
{
    FieldManager *fm = this->domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf;
    int err;
    if ( ( tf = fm->giveField( FT_Temperature ) ) ) {
        // temperature field registered
        FloatArray gcoords, answer;
        static_cast<StructuralElement *>( gp->giveElement() )->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
        if ( ( err = tf->evaluateAt( answer, gcoords, VM_Total, tStep ) ) ) {
            OOFEM_ERROR( "tf->evaluateAt failed, element %d, error code %d", gp->giveElement()->giveNumber(), err );
        }
        return answer.at( 1 );
    }
    return 0.;
}

double IsotropicLinearElasticMaterialTemperature::eval_E( GaussPoint *gp, TimeStep *tStep ) const
{
    const std::map<std::string, FunctionArgument> MapTemp{
        { "te", giveTemperature( gp, tStep ) },
    };

    return E.eval( MapTemp, this->giveDomain(), gp, giveTemperature( gp, tStep ) );
}


double IsotropicLinearElasticMaterialTemperature::eval_nu( GaussPoint *gp, TimeStep *tStep ) const
{
    const std::map<std::string, FunctionArgument> MapTemp{
        { "te", giveTemperature( gp, tStep ) },
    };

    return nu.eval( MapTemp, this->giveDomain(), gp, giveTemperature( gp, tStep ) );
}
} // end namespace oofem
