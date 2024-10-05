//
// Created by babaid on 05.10.24.
//

#include "openmm/common/PuremdInterface.h"
#include "openmm/OpenMMException.h"

#include "spuremd.h"

using namespace OpenMM;
PuremdInterface::PuremdInterface( const std::vector<double> & sim_box_info, const std::string &ffield_filename, const std::string &control_filename): firstCall(true), sim_box_info(sim_box_info), ffield_filename(ffield_filename), control_filename(control_filename) {}

void PuremdInterface::getReaxffPuremdForces(int num_qm_atoms, std::string &qm_symbols, const std::vector<double> & qm_pos,
                                            int num_mm_atoms, const std::string &mm_symbols, const std::vector<double> & mm_pos_q,
                                            std::vector<double>& qm_forces, std::vector<double>& mm_forces, std::vector<double> qm_q, double& totalEnergy) {
  qm_forces.resize(num_qm_atoms);
  mm_forces.resize(num_mm_atoms);
  qm_q.resize(num_qm_atoms);


  if(firstCall)
  {
    if (!control_filename.empty()) {
      handlePuremd = setup_qmmm(
          num_qm_atoms, qm_symbols.data(), qm_pos.data(), num_mm_atoms,
          mm_symbols.data(), mm_pos_q.data(), sim_box_info.data(),
          ffield_filename.data(), control_filename.data());
    }
    else
    {
      handlePuremd = setup_qmmm(
          num_qm_atoms, qm_symbols.data(), qm_pos.data(), num_mm_atoms,
          mm_symbols.data(), mm_pos_q.data(), sim_box_info.data(),
          ffield_filename.data(), NULL);

      retPuremd = simulate(handlePuremd);
      if (0 != retPuremd) throw OpenMMException("Error at PuReMD simulation.");


    }
    firstCall = false;
  }
  else
  {
      retPuremd = reset_qmmm(handlePuremd, num_qm_atoms, qm_symbols.data(), qm_pos.data(),
                     num_mm_atoms, mm_symbols.data(), mm_pos_q.data(),
                     sim_box_info.data(),
                     NULL, NULL);
      if(0 != retPuremd) throw OpenMMException("Issue with PuReMD function reset_qmmm.");
  }

  retPuremd = get_atom_forces_qmmm(handlePuremd, qm_forces.data(), mm_forces.data());
  retPuremd = get_atom_charges_qmmm(handlePuremd, qm_q.data(), NULL);
  retPuremd = get_system_info(handlePuremd, NULL, NULL, &totalEnergy, NULL, NULL, NULL);


}

