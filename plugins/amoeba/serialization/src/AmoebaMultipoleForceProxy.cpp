/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/serialization/AmoebaMultipoleForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaMultipoleForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaMultipoleForceProxy::AmoebaMultipoleForceProxy() : SerializationProxy("AmoebaMultipoleForce") {
}

static void getCovalentTypes( std::vector<std::string>& covalentTypes ){

    covalentTypes.push_back( "Covalent12" );
    covalentTypes.push_back( "Covalent13" );
    covalentTypes.push_back( "Covalent14" );
    covalentTypes.push_back( "Covalent15" );

    covalentTypes.push_back( "PolarizationCovalent11" );
    covalentTypes.push_back( "PolarizationCovalent12" );
    covalentTypes.push_back( "PolarizationCovalent13" );
    covalentTypes.push_back( "PolarizationCovalent14" );
}

static void addCovalentMap( SerializationNode& particleExclusions, int particleIndex, std::string mapName, std::vector< int > covalentMap ){
    SerializationNode& map   = particleExclusions.createChildNode(mapName);
    for (unsigned int ii = 0; ii < covalentMap.size(); ii++) {
        map.createChildNode("Cv").setIntProperty( "v", covalentMap[ii] );
    }
}

void loadCovalentMap( const SerializationNode& map, std::vector< int >& covalentMap ){
    for (unsigned int ii = 0; ii < map.getChildren().size(); ii++) {
        covalentMap.push_back( map.getChildren()[ii].getIntProperty( "v" ) );
    }
}

void AmoebaMultipoleForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const AmoebaMultipoleForce& force = *reinterpret_cast<const AmoebaMultipoleForce*>(object);

    node.setIntProperty("nonbondedMethod",                  force.getNonbondedMethod());
    node.setIntProperty("polarizationType",                 force.getPolarizationType());
    //node.setIntProperty("pmeBSplineOrder",                  force.getPmeBSplineOrder());
    //node.setIntProperty("mutualInducedIterationMethod",     force.getMutualInducedIterationMethod());
    node.setIntProperty("mutualInducedMaxIterations",       force.getMutualInducedMaxIterations());

    node.setDoubleProperty("cutoffDistance",                force.getCutoffDistance());
    node.setDoubleProperty("aEwald",                        force.getAEwald());
    node.setDoubleProperty("mutualInducedTargetEpsilon",    force.getMutualInducedTargetEpsilon());
    //node.setDoubleProperty("electricConstant",              force.getElectricConstant());
    node.setDoubleProperty("ewaldErrorTolerance",           force.getEwaldErrorTolerance());

    std::vector<int> gridDimensions;
    force.getPmeGridDimensions( gridDimensions );
    SerializationNode& gridDimensionsNode  = node.createChildNode("MultipoleParticleGridDimension");
    gridDimensionsNode.setIntProperty( "d0", gridDimensions[0] ).setIntProperty( "d1", gridDimensions[1] ).setIntProperty( "d2", gridDimensions[2] ); 

    std::vector<std::string> covalentTypes;
    getCovalentTypes( covalentTypes );

    SerializationNode& particles = node.createChildNode("MultipoleParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumMultipoles()); ii++) {

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, thole, dampingFactor, polarity;

        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;

        force.getMultipoleParameters( ii, charge, molecularDipole, molecularQuadrupole,
                                      axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity );

        SerializationNode& particle    = particles.createChildNode("Particle");
        particle.setIntProperty("axisType", axisType).setIntProperty("multipoleAtomZ", multipoleAtomZ).setIntProperty("multipoleAtomX", multipoleAtomX).setIntProperty("multipoleAtomY", multipoleAtomY);
        particle.setDoubleProperty("charge", charge).setDoubleProperty("thole", thole).setDoubleProperty("damp", dampingFactor).setDoubleProperty("polarity", polarity);

        SerializationNode& dipole      = particle.createChildNode("Dipole");
        dipole.setDoubleProperty( "d0", molecularDipole[0] ).setDoubleProperty( "d1", molecularDipole[1] ).setDoubleProperty( "d2", molecularDipole[2] );

        SerializationNode& quadrupole  = particle.createChildNode("Quadrupole");
        quadrupole.setDoubleProperty( "q0", molecularQuadrupole[0] ).setDoubleProperty( "q1", molecularQuadrupole[1] ).setDoubleProperty( "q2", molecularQuadrupole[2] );
        quadrupole.setDoubleProperty( "q3", molecularQuadrupole[3] ).setDoubleProperty( "q4", molecularQuadrupole[4] ).setDoubleProperty( "q5", molecularQuadrupole[5] );
        quadrupole.setDoubleProperty( "q6", molecularQuadrupole[6] ).setDoubleProperty( "q7", molecularQuadrupole[7] ).setDoubleProperty( "q8", molecularQuadrupole[8] );

        for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
            std::vector< int > covalentMap;
            force.getCovalentMap(ii, static_cast<AmoebaMultipoleForce::CovalentType>(jj), covalentMap );
            addCovalentMap( particle, ii, covalentTypes[jj], covalentMap );
        }
    }
}

void* AmoebaMultipoleForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") > 2)
        throw OpenMMException("Unsupported version number");
    AmoebaMultipoleForce* force = new AmoebaMultipoleForce();

    try {

        force->setNonbondedMethod( static_cast<AmoebaMultipoleForce::NonbondedMethod>(node.getIntProperty( "nonbondedMethod" )) );
        if( node.getIntProperty("version") == 2 ){
            force->setPolarizationType( static_cast<AmoebaMultipoleForce::PolarizationType>(node.getIntProperty( "polarizationType" )) );
        }
        //force->setPmeBSplineOrder( node.getIntProperty( "pmeBSplineOrder" ) );
        //force->setMutualInducedIterationMethod( static_cast<AmoebaMultipoleForce::MutualInducedIterationMethod>(node.getIntProperty( "mutualInducedIterationMethod" ) ) );
        force->setMutualInducedMaxIterations( node.getIntProperty( "mutualInducedMaxIterations" ) );

        force->setCutoffDistance( node.getDoubleProperty( "cutoffDistance" ) );
        force->setAEwald( node.getDoubleProperty( "aEwald" ) );
        force->setMutualInducedTargetEpsilon( node.getDoubleProperty( "mutualInducedTargetEpsilon" ) );
        //force->setElectricConstant( node.getDoubleProperty( "electricConstant" ) );
        force->setEwaldErrorTolerance( node.getDoubleProperty( "ewaldErrorTolerance" ) );

        std::vector<int> gridDimensions;
        const SerializationNode& gridDimensionsNode  = node.getChildNode("MultipoleParticleGridDimension");
        gridDimensions.push_back( gridDimensionsNode.getIntProperty( "d0" ));
        gridDimensions.push_back( gridDimensionsNode.getIntProperty( "d1" ));
        gridDimensions.push_back( gridDimensionsNode.getIntProperty( "d2" ));
        force->setPmeGridDimensions( gridDimensions );
    
        std::vector<std::string> covalentTypes;
        getCovalentTypes( covalentTypes );

        const SerializationNode& particles = node.getChildNode("MultipoleParticles");
        for ( unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {

            const SerializationNode& particle = particles.getChildren()[ii];

            std::vector<double> molecularDipole;
            const SerializationNode& dipole = particle.getChildNode("Dipole");
            molecularDipole.push_back( dipole.getDoubleProperty( "d0" ) );
            molecularDipole.push_back( dipole.getDoubleProperty( "d1" ) );
            molecularDipole.push_back( dipole.getDoubleProperty( "d2" ) );

            std::vector<double> molecularQuadrupole;
            const SerializationNode& quadrupole = particle.getChildNode("Quadrupole");
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q0" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q1" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q2" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q3" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q4" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q5" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q6" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q7" ) );
            molecularQuadrupole.push_back( quadrupole.getDoubleProperty( "q8" ) );

            force->addMultipole( particle.getDoubleProperty("charge"), molecularDipole, molecularQuadrupole,
                                particle.getIntProperty("axisType"),
                                particle.getIntProperty("multipoleAtomZ"),
                                particle.getIntProperty("multipoleAtomX"),
                                particle.getIntProperty("multipoleAtomY"),
                                particle.getDoubleProperty("thole"),
                                particle.getDoubleProperty("damp"), particle.getDoubleProperty("polarity"));

            // covalent maps 

            for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
                std::vector< int > covalentMap;
                loadCovalentMap( particle.getChildNode(covalentTypes[jj]), covalentMap );
                force->setCovalentMap( ii, static_cast<AmoebaMultipoleForce::CovalentType>(jj), covalentMap );
            }
        }

    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}
