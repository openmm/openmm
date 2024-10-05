//
// Created by babaid on 05.10.24.
//


#ifndef OPENMM_PUREMDINTERFACE_H
#define OPENMM_PUREMDINTERFACE_H
#define QMMM

#include<vector>
#include<string>

class PuremdInterface
{
private:
  bool firstCall;
  void* handlePuremd;
  int retPuremd;
  const std::vector<double>  sim_box_info;
  const std::string ffield_filename;
  const std::string control_filename;
public:
  PuremdInterface( const std::vector<double> & sim_box_info, const std::string &ffield_filename, const std::string &control_filename);
  void getReaxffPuremdForces(int num_qm_atoms, std::string &qm_symbols, const std::vector<double> & qm_pos,
                             int num_mm_atoms, const std::string &mm_symbols, const std::vector<double> & mm_pos_q,
                             std::vector<double>& qm_forces, std::vector<double>& mm_forces, std::vector<double> qm_q, double& totalEnergy);
};

#endif // OPENMM_PUREMDINTERFACE_H
