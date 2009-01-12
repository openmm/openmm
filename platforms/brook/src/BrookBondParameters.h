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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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
