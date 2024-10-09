//
// Created by babaid on 05.10.24.
//

#include "openmm/common/PuremdInterface.h"
#include "openmm/OpenMMException.h"

#include "spuremd.h"

using namespace OpenMM;
PuremdInterface::PuremdInterface(): firstCall(true) {}


void PuremdInterface::setInputFileNames(const std::string &ffieldFilename, const std::string &controlFilename) {

  ffield_filename = ffieldFilename;
  control_filename = controlFilename;
}

void PuremdInterface::getReaxffPuremdForces(int num_qm_atoms,  const std::vector<char> &qm_symbols, const std::vector<double> & qm_pos,
                                            int num_mm_atoms, const  std::vector<char> &mm_symbols, const std::vector<double> & mm_pos_q,
                                            const std::vector<double> & sim_box_info,
                                            std::vector<double>& new_qm_pos, std::vector<double>& new_mm_pos,
                                            std::vector<double>& qm_forces, std::vector<double>& mm_forces, std::vector<double> qm_q, double& totalEnergy) {

  if(firstCall)
  {
    if (!control_filename.empty()) {
      handlePuremd = setup_qmmm(
          num_qm_atoms, qm_symbols.data(), qm_pos.data(), num_mm_atoms,
          mm_symbols.data(), mm_pos_q.data(), sim_box_info.data(),
          ffield_filename.c_str(), control_filename.c_str());
    }
    else
    {
      handlePuremd = setup_qmmm(
          num_qm_atoms, qm_symbols.data(), qm_pos.data(), num_mm_atoms,
          mm_symbols.data(), mm_pos_q.data(), sim_box_info.data(),
          ffield_filename.c_str(), NULL);





    }
    firstCall = false;
  }
  else
  {
      retPuremd = reset_qmmm(handlePuremd, num_qm_atoms, qm_symbols.data(),  qm_pos.data(),
                     num_mm_atoms, mm_symbols.data(), mm_pos_q.data(),
                     sim_box_info.data(),
                     NULL, NULL);
      if(0 != retPuremd) throw OpenMMException("Issue with PuReMD function reset_qmmm.");
  }

  retPuremd = simulate(handlePuremd);
  if (0 != retPuremd) throw OpenMMException("Error at PuReMD simulation.");
  retPuremd = get_atom_forces_qmmm(handlePuremd, qm_forces.data(), mm_forces.data());
  retPuremd = get_atom_charges_qmmm(handlePuremd, qm_q.data(), NULL);
  retPuremd = get_system_info(handlePuremd, NULL, NULL, &totalEnergy, NULL, NULL, NULL);
  retPuremd = get_atom_positions_qmmm(handlePuremd, new_qm_pos.data(), new_mm_pos.data());
}

