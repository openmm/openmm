#ifndef OPENMM_BROOK_UPDATE_TIME_KERNEL_H_
#define OPENMM_BROOK_UPDATE_TIME_KERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston, Peter Eastman                      *
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

#include "openmm/kernels.h"
#include "OpenMMBrookInterface.h"

namespace OpenMM {

/**
 * This kernel initializes the forces
 */
class BrookUpdateTimeKernel : public UpdateTimeKernel {

   public:

      BrookUpdateTimeKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface );

      ~BrookUpdateTimeKernel();

      /**
       * Initialize the kernel
       *
       * @param system     the System this kernel will be applied to
       */

      void initialize( const System& system );

      /**
       * Get the current time (in picoseconds).
       *
       * @param context    the context in which to execute this kernel
       */

      double getTime(const ContextImpl& context) const;

      /**
       * Set the current time (in picoseconds).
       *
       * @param context    the context in which to execute this kernel
       * @param time       the time
       */

      void setTime(ContextImpl& context, double time);


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

   private:

      // log file reference

      FILE* _log;

      // interface

      OpenMMBrookInterface& _openMMBrookInterface;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_UPDATE_TIME_KERNEL_H_ */

