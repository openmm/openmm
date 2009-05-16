#ifndef OPENMM_BROOK_BOND_PARAMETERS_H_
#define OPENMM_BROOK_BOND_PARAMETERS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include <vector>

namespace OpenMM {

/**
 * Container for bond parameters
 */
class BrookBondParameters {

   public:
  
     // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1; 

      /** 
       * BrookBondParameters constructor
       * 
       * @param bondName                  bond name
       * @param numberOfParticlesInBond   no. of particles in each bond
       * @param numberOfParametersInBond  no. of parameters in each bond
       * @param numberOfBonds             no. of bonds
       * @param log                       optional log reference
       *
       */
      
      BrookBondParameters( std::string bondName, int numberOfParticlesInBond, int numberOfParametersInBond, int numberOfBonds, FILE* log );
  
      ~BrookBondParameters();
  
      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
      /** 
       * Set bond info
       * 
       * @param bondIndex       index of bond
       * @param particleIndices array of particle indices
       * @param bondParameters  array of bond parameters
       *
       * @return  DefaultReturnValue
       *
       * @throw OpenMMException exeception if bond index is invalid
       *
       */
      
      int setBond( int bondIndex, int* particleIndices, double* bondParameters );
      
      /** 
       * Get bond name
       * 
       * @return  bond name
       *
       */
      
      std::string getBondName( void ) const;
      
      /** 
       * Set bond name
       * 
       * @param bondName       bond name
       *
       * @return  DefaultReturnValue
       *
       */
      
      //int setBondName( std::string bondName );
      
      /** 
       * Get NumberOfParticlesInBond
       * 
       * @return NumberOfParticlesInBond
       *
       */
      
      int getNumberOfParticlesInBond( void ) const;
      
      /** 
       * Get NumberOfParametersInBond
       * 
       * @return NumberOfParametersInBond
       *
       */
      
      int getNumberOfParametersInBond( void ) const;
      
      /** 
       * Get NumberOfBonds
       * 
       * @return NumberOfBonds
       *
       */
      
      int getNumberOfBonds( void ) const;
      
      /** 
       * Get particle indices
       * 
       * @return particle indices
       *
       */
      
      const std::vector<std::vector<int>>& getParticleIndices( void ) const;
      
      /** 
       * Get parameters
       * 
       * @return parameters
       *
       */
      
      const std::vector<std::vector<double>>& getBondParameters( void ) const;
      
      /*  
       * Get contents of object
       *
       *
       * @param level   level of dump
       *
       * @return string containing contents
       *
       * */
    
      std::string getContentsString( int level = 0 ) const;

   private:
   
      // log file reference

      FILE* _log;

      // bond name

      std::string _bondName;

      // number of bonds

      int _numberOfBonds;
      int _numberOfParticlesInBond;
      int _numberOfParametersInBond;
   
      // particle indices and parameters

      std::vector<std::vector<int>>    _particleIndices;
      std::vector<std::vector<double>> _bondParameters;

      /* 
       * Get contents of object
       *
       * @param tab         tab
       * @param description description
       * @param value       value
       *
       * @return string containing contents
       *
       * */
      
      std::string _getLine( const std::string& tab, const std::string& description,
                            const std::string& value ) const;
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_BOND_PARAMETERS_H_ */
