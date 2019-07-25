
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __ObcParameters_H__
#define __ObcParameters_H__

#include "RealVec.h"
#include <vector>

namespace OpenMM {

class ObcParameters {

   public:

       // OBC types

       enum ObcType { ObcTypeI, ObcTypeII };

   private:

      // OBC constants & parameters
   
      int _numberOfAtoms;

      RealOpenMM _solventDielectric;
      RealOpenMM _soluteDielectric;
      RealOpenMM _electricConstant;
      RealOpenMM _probeRadius;
      RealOpenMM _pi4Asolv;

      RealOpenMM _dielectricOffset;
      RealOpenMM _alphaObc;
      RealOpenMM _betaObc;
      RealOpenMM _gammaObc;
      ObcType _obcType;

      // scaled radius factors (S_kk in HCT paper)

      std::vector<RealOpenMM> _atomicRadii;
      std::vector<RealOpenMM> _scaledRadiusFactors;

      // cutoff and periodic boundary conditions
      
      bool _cutoff;
      bool _periodic;
      OpenMM::RealVec _periodicBoxVectors[3];
      RealOpenMM _cutoffDistance;

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric 
      
         @param dielectricOffset         solvent dielectric

         --------------------------------------------------------------------------------------- */
      
      void setDielectricOffset(RealOpenMM dielectricOffset);

   public:

      /**---------------------------------------------------------------------------------------
      
         ObcParameters constructor 
      
         @param numberOfAtoms       number of atoms
      
         --------------------------------------------------------------------------------------- */
      
       ObcParameters(int numberOfAtoms, ObcParameters::ObcType obcType = ObcTypeII);

      /**---------------------------------------------------------------------------------------
      
         ObcParameters destructor 
      
         --------------------------------------------------------------------------------------- */
      
       ~ObcParameters();

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */

      int getNumberOfAtoms() const;

      /**---------------------------------------------------------------------------------------
      
         Get electric constant
      
         @return electric constant
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getElectricConstant() const;

      /**---------------------------------------------------------------------------------------
      
         Get probe radius (Simbios) 
      
         @return probeRadius
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getProbeRadius() const;

      /**---------------------------------------------------------------------------------------
      
         Set probe radius (Simbios) 
      
         @param probeRadius   probe radius
      
         --------------------------------------------------------------------------------------- */

      void setProbeRadius(RealOpenMM probeRadius);

      /**---------------------------------------------------------------------------------------
      
         Get pi4Asolv:  used in ACE approximation for nonpolar term  
            ((RealOpenMM) M_PI)*4.0f*0.0049f*1000.0f; (Simbios) 
      
         @return pi4Asolv
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getPi4Asolv() const;

      /**---------------------------------------------------------------------------------------
      
         Set pi4Asolv
      
         --------------------------------------------------------------------------------------- */

      void setPi4Asolv(RealOpenMM pi4Asolv);

      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric
      
         @return solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSolventDielectric() const;

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric
      
         @param solventDielectric solvent dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSolventDielectric(RealOpenMM solventDielectric);

      /**---------------------------------------------------------------------------------------
      
         Get solute dielectric
      
         @return soluteDielectric
      
         --------------------------------------------------------------------------------------- */

      RealOpenMM getSoluteDielectric() const;

      /**---------------------------------------------------------------------------------------
      
         Set solute dielectric
      
         @param soluteDielectric solute dielectric
      
         --------------------------------------------------------------------------------------- */

      void setSoluteDielectric(RealOpenMM soluteDielectric);

      /**---------------------------------------------------------------------------------------
      
         Get OBC type
      
         @return OBC type
      
         --------------------------------------------------------------------------------------- */
      
      ObcParameters::ObcType getObcType() const;
      
      /**---------------------------------------------------------------------------------------
      
         Set OBC type specific parameters
      
         @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)
      
         --------------------------------------------------------------------------------------- */
      
      void setObcTypeParameters(ObcParameters::ObcType obcType);
      
      /**---------------------------------------------------------------------------------------
      
         Get alpha OBC (Eqs. 6 & 7) in Proteins paper
      
         @return alphaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getAlphaObc() const;
      
      /**---------------------------------------------------------------------------------------
      
         Get beta OBC (Eqs. 6 & 7) in Proteins paper
      
         @return betaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getBetaObc() const;
      
      /**---------------------------------------------------------------------------------------
      
         Get gamma OBC (Eqs. 6 & 7) in Proteins paper
      
         @return gammaObc
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getGammaObc() const;
      
      /**---------------------------------------------------------------------------------------
      
         Get solvent dielectric 
      
         @return dielectricOffset dielectric offset
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getDielectricOffset() const;

      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const std::vector<RealOpenMM>& getScaledRadiusFactors() const;
        
      /**---------------------------------------------------------------------------------------
      
         Return OBC scale factors
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      void setScaledRadiusFactors(const std::vector<RealOpenMM>& scaledRadiusFactors);
        
      /**---------------------------------------------------------------------------------------
      
         Get AtomicRadii array w/ dielectric offset applied
      
         @return array of atom volumes
      
         --------------------------------------------------------------------------------------- */

      const std::vector<RealOpenMM>& getAtomicRadii() const;

      /**---------------------------------------------------------------------------------------
      
         Set AtomicRadii array
      
         @param atomicRadii vector of atomic radii
      
         --------------------------------------------------------------------------------------- */

      void setAtomicRadii(const std::vector<RealOpenMM>& atomicRadii);


      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(RealOpenMM distance);

      /**---------------------------------------------------------------------------------------

         Get whether to use a cutoff.

         --------------------------------------------------------------------------------------- */

      bool getUseCutoff() const;

      /**---------------------------------------------------------------------------------------

         Get the cutoff distance.

         --------------------------------------------------------------------------------------- */

      RealOpenMM getCutoffDistance() const;

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param vectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(OpenMM::RealVec* vectors);

      /**---------------------------------------------------------------------------------------

         Get whether to use periodic boundary conditions.

         --------------------------------------------------------------------------------------- */

      bool getPeriodic();

      /**---------------------------------------------------------------------------------------

         Get the periodic box vectors

         --------------------------------------------------------------------------------------- */

      const OpenMM::RealVec* getPeriodicBox();

};

} // namespace OpenMM

#endif // __ObcParameters_H__
