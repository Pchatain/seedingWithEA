#include "ActsExamples/Io/Root/RootLdmxSimHitReader.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

// The Ldmx space point to save into the eventStore
#include "Acts/LdmxEDM/LdmxSpacePoint.hpp"

#include <iostream>

#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>

ActsExamples::RootLdmxSimHitReader::RootLdmxSimHitReader(
    ActsExamples::RootLdmxSimHitReader::Config cfg, Acts::Logging::Level lvl)
    : m_cfg(std::move(cfg)),
      m_events(0),
      fChain(nullptr),
      m_logger(Acts::getDefaultLogger("RootLdmxSimHitReader", lvl)) {
  ACTS_DEBUG("Constructor of RootLdmxSimHitReader");
  ACTS_DEBUG("The tree name" << m_cfg.treeName);

  fChain = new TChain(m_cfg.treeName.c_str());

  ACTS_DEBUG("fChain = " << fChain);

  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("SimParticles_v12", &SimParticles_v12_,
                           &b_SimParticles_v12_);
  fChain->SetBranchAddress("SimParticles_v12.first", SimParticles_v12_first,
                           &b_SimParticles_v12_first);
  fChain->SetBranchAddress("SimParticles_v12.second.energy_",
                           SimParticles_v12_second_energy_,
                           &b_SimParticles_v12_second_energy_);
  fChain->SetBranchAddress("SimParticles_v12.second.pdgID_",
                           SimParticles_v12_second_pdgID_,
                           &b_SimParticles_v12_second_pdgID_);
  fChain->SetBranchAddress("SimParticles_v12.second.genStatus_",
                           SimParticles_v12_second_genStatus_,
                           &b_SimParticles_v12_second_genStatus_);
  fChain->SetBranchAddress("SimParticles_v12.second.time_",
                           SimParticles_v12_second_time_,
                           &b_SimParticles_v12_second_time_);
  fChain->SetBranchAddress("SimParticles_v12.second.x_",
                           SimParticles_v12_second_x_,
                           &b_SimParticles_v12_second_x_);
  fChain->SetBranchAddress("SimParticles_v12.second.y_",
                           SimParticles_v12_second_y_,
                           &b_SimParticles_v12_second_y_);
  fChain->SetBranchAddress("SimParticles_v12.second.z_",
                           SimParticles_v12_second_z_,
                           &b_SimParticles_v12_second_z_);
  fChain->SetBranchAddress("SimParticles_v12.second.endX_",
                           SimParticles_v12_second_endX_,
                           &b_SimParticles_v12_second_endX_);
  fChain->SetBranchAddress("SimParticles_v12.second.endY_",
                           SimParticles_v12_second_endY_,
                           &b_SimParticles_v12_second_endY_);
  fChain->SetBranchAddress("SimParticles_v12.second.endZ_",
                           SimParticles_v12_second_endZ_,
                           &b_SimParticles_v12_second_endZ_);
  fChain->SetBranchAddress("SimParticles_v12.second.px_",
                           SimParticles_v12_second_px_,
                           &b_SimParticles_v12_second_px_);
  fChain->SetBranchAddress("SimParticles_v12.second.py_",
                           SimParticles_v12_second_py_,
                           &b_SimParticles_v12_second_py_);
  fChain->SetBranchAddress("SimParticles_v12.second.pz_",
                           SimParticles_v12_second_pz_,
                           &b_SimParticles_v12_second_pz_);
  fChain->SetBranchAddress("SimParticles_v12.second.endpx_",
                           SimParticles_v12_second_endpx_,
                           &b_SimParticles_v12_second_endpx_);
  fChain->SetBranchAddress("SimParticles_v12.second.endpy_",
                           SimParticles_v12_second_endpy_,
                           &b_SimParticles_v12_second_endpy_);
  fChain->SetBranchAddress("SimParticles_v12.second.endpz_",
                           SimParticles_v12_second_endpz_,
                           &b_SimParticles_v12_second_endpz_);
  fChain->SetBranchAddress("SimParticles_v12.second.mass_",
                           SimParticles_v12_second_mass_,
                           &b_SimParticles_v12_second_mass_);
  fChain->SetBranchAddress("SimParticles_v12.second.charge_",
                           SimParticles_v12_second_charge_,
                           &b_SimParticles_v12_second_charge_);
  fChain->SetBranchAddress("SimParticles_v12.second.daughters_",
                           SimParticles_v12_second_daughters_,
                           &b_SimParticles_v12_second_daughters_);
  fChain->SetBranchAddress("SimParticles_v12.second.parents_",
                           SimParticles_v12_second_parents_,
                           &b_SimParticles_v12_second_parents_);
  fChain->SetBranchAddress("SimParticles_v12.second.processType_",
                           SimParticles_v12_second_processType_,
                           &b_SimParticles_v12_second_processType_);
  fChain->SetBranchAddress("SimParticles_v12.second.vertexVolume_",
                           SimParticles_v12_second_vertexVolume_,
                           &b_SimParticles_v12_second_vertexVolume_);
  fChain->SetBranchAddress("TaggerSimHits_v12", &TaggerSimHits_v12_,
                           &b_TaggerSimHits_v12_);
  fChain->SetBranchAddress("TaggerSimHits_v12.id_", TaggerSimHits_v12_id_,
                           &b_TaggerSimHits_v12_id_);
  fChain->SetBranchAddress("TaggerSimHits_v12.layerID_",
                           TaggerSimHits_v12_layerID_,
                           &b_TaggerSimHits_v12_layerID_);
  fChain->SetBranchAddress("TaggerSimHits_v12.moduleID_",
                           TaggerSimHits_v12_moduleID_,
                           &b_TaggerSimHits_v12_moduleID_);
  fChain->SetBranchAddress("TaggerSimHits_v12.edep_", TaggerSimHits_v12_edep_,
                           &b_TaggerSimHits_v12_edep_);
  fChain->SetBranchAddress("TaggerSimHits_v12.time_", TaggerSimHits_v12_time_,
                           &b_TaggerSimHits_v12_time_);
  fChain->SetBranchAddress("TaggerSimHits_v12.px_", TaggerSimHits_v12_px_,
                           &b_TaggerSimHits_v12_px_);
  fChain->SetBranchAddress("TaggerSimHits_v12.py_", TaggerSimHits_v12_py_,
                           &b_TaggerSimHits_v12_py_);
  fChain->SetBranchAddress("TaggerSimHits_v12.pz_", TaggerSimHits_v12_pz_,
                           &b_TaggerSimHits_v12_pz_);
  fChain->SetBranchAddress("TaggerSimHits_v12.energy_",
                           TaggerSimHits_v12_energy_,
                           &b_TaggerSimHits_v12_energy_);
  fChain->SetBranchAddress("TaggerSimHits_v12.x_", TaggerSimHits_v12_x_,
                           &b_TaggerSimHits_v12_x_);
  fChain->SetBranchAddress("TaggerSimHits_v12.y_", TaggerSimHits_v12_y_,
                           &b_TaggerSimHits_v12_y_);
  fChain->SetBranchAddress("TaggerSimHits_v12.z_", TaggerSimHits_v12_z_,
                           &b_TaggerSimHits_v12_z_);
  fChain->SetBranchAddress("TaggerSimHits_v12.pathLength_",
                           TaggerSimHits_v12_pathLength_,
                           &b_TaggerSimHits_v12_pathLength_);
  fChain->SetBranchAddress(
      "TaggerSimHits_v12.trackID_", TaggerSimHits_v12_trackID_,
      &b_TaggerSimHits_v12_trackID_);  // is this the particle info for this
                                       // hit?
  fChain->SetBranchAddress("TaggerSimHits_v12.pdgID_", TaggerSimHits_v12_pdgID_,
                           &b_TaggerSimHits_v12_pdgID_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12", &TriggerPadUpSimHits_v12_,
                           &b_TriggerPadUpSimHits_v12_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.id_",
                           TriggerPadUpSimHits_v12_id_,
                           &b_TriggerPadUpSimHits_v12_id_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.edep_",
                           TriggerPadUpSimHits_v12_edep_,
                           &b_TriggerPadUpSimHits_v12_edep_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.x_",
                           TriggerPadUpSimHits_v12_x_,
                           &b_TriggerPadUpSimHits_v12_x_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.y_",
                           TriggerPadUpSimHits_v12_y_,
                           &b_TriggerPadUpSimHits_v12_y_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.z_",
                           TriggerPadUpSimHits_v12_z_,
                           &b_TriggerPadUpSimHits_v12_z_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.time_",
                           TriggerPadUpSimHits_v12_time_,
                           &b_TriggerPadUpSimHits_v12_time_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.trackIDContribs_",
                           TriggerPadUpSimHits_v12_trackIDContribs_,
                           &b_TriggerPadUpSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.incidentIDContribs_",
                           TriggerPadUpSimHits_v12_incidentIDContribs_,
                           &b_TriggerPadUpSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.pdgCodeContribs_",
                           TriggerPadUpSimHits_v12_pdgCodeContribs_,
                           &b_TriggerPadUpSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.edepContribs_",
                           TriggerPadUpSimHits_v12_edepContribs_,
                           &b_TriggerPadUpSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.timeContribs_",
                           TriggerPadUpSimHits_v12_timeContribs_,
                           &b_TriggerPadUpSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("TriggerPadUpSimHits_v12.nContribs_",
                           TriggerPadUpSimHits_v12_nContribs_,
                           &b_TriggerPadUpSimHits_v12_nContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12",
                           &TriggerPadTaggerSimHits_v12_,
                           &b_TriggerPadTaggerSimHits_v12_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.id_",
                           TriggerPadTaggerSimHits_v12_id_,
                           &b_TriggerPadTaggerSimHits_v12_id_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.edep_",
                           TriggerPadTaggerSimHits_v12_edep_,
                           &b_TriggerPadTaggerSimHits_v12_edep_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.x_",
                           TriggerPadTaggerSimHits_v12_x_,
                           &b_TriggerPadTaggerSimHits_v12_x_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.y_",
                           TriggerPadTaggerSimHits_v12_y_,
                           &b_TriggerPadTaggerSimHits_v12_y_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.z_",
                           TriggerPadTaggerSimHits_v12_z_,
                           &b_TriggerPadTaggerSimHits_v12_z_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.time_",
                           TriggerPadTaggerSimHits_v12_time_,
                           &b_TriggerPadTaggerSimHits_v12_time_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.trackIDContribs_",
                           TriggerPadTaggerSimHits_v12_trackIDContribs_,
                           &b_TriggerPadTaggerSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.incidentIDContribs_",
                           TriggerPadTaggerSimHits_v12_incidentIDContribs_,
                           &b_TriggerPadTaggerSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.pdgCodeContribs_",
                           TriggerPadTaggerSimHits_v12_pdgCodeContribs_,
                           &b_TriggerPadTaggerSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.edepContribs_",
                           TriggerPadTaggerSimHits_v12_edepContribs_,
                           &b_TriggerPadTaggerSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.timeContribs_",
                           TriggerPadTaggerSimHits_v12_timeContribs_,
                           &b_TriggerPadTaggerSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("TriggerPadTaggerSimHits_v12.nContribs_",
                           TriggerPadTaggerSimHits_v12_nContribs_,
                           &b_TriggerPadTaggerSimHits_v12_nContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12", &TargetSimHits_v12_,
                           &b_TargetSimHits_v12_);
  fChain->SetBranchAddress("TargetSimHits_v12.id_", TargetSimHits_v12_id_,
                           &b_TargetSimHits_v12_id_);
  fChain->SetBranchAddress("TargetSimHits_v12.edep_", TargetSimHits_v12_edep_,
                           &b_TargetSimHits_v12_edep_);
  fChain->SetBranchAddress("TargetSimHits_v12.x_", TargetSimHits_v12_x_,
                           &b_TargetSimHits_v12_x_);
  fChain->SetBranchAddress("TargetSimHits_v12.y_", TargetSimHits_v12_y_,
                           &b_TargetSimHits_v12_y_);
  fChain->SetBranchAddress("TargetSimHits_v12.z_", TargetSimHits_v12_z_,
                           &b_TargetSimHits_v12_z_);
  fChain->SetBranchAddress("TargetSimHits_v12.time_", TargetSimHits_v12_time_,
                           &b_TargetSimHits_v12_time_);
  fChain->SetBranchAddress("TargetSimHits_v12.trackIDContribs_",
                           TargetSimHits_v12_trackIDContribs_,
                           &b_TargetSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12.incidentIDContribs_",
                           TargetSimHits_v12_incidentIDContribs_,
                           &b_TargetSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12.pdgCodeContribs_",
                           TargetSimHits_v12_pdgCodeContribs_,
                           &b_TargetSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12.edepContribs_",
                           TargetSimHits_v12_edepContribs_,
                           &b_TargetSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12.timeContribs_",
                           TargetSimHits_v12_timeContribs_,
                           &b_TargetSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("TargetSimHits_v12.nContribs_",
                           TargetSimHits_v12_nContribs_,
                           &b_TargetSimHits_v12_nContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12",
                           &TriggerPadDownSimHits_v12_,
                           &b_TriggerPadDownSimHits_v12_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.id_",
                           TriggerPadDownSimHits_v12_id_,
                           &b_TriggerPadDownSimHits_v12_id_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.edep_",
                           TriggerPadDownSimHits_v12_edep_,
                           &b_TriggerPadDownSimHits_v12_edep_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.x_",
                           TriggerPadDownSimHits_v12_x_,
                           &b_TriggerPadDownSimHits_v12_x_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.y_",
                           TriggerPadDownSimHits_v12_y_,
                           &b_TriggerPadDownSimHits_v12_y_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.z_",
                           TriggerPadDownSimHits_v12_z_,
                           &b_TriggerPadDownSimHits_v12_z_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.time_",
                           TriggerPadDownSimHits_v12_time_,
                           &b_TriggerPadDownSimHits_v12_time_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.trackIDContribs_",
                           TriggerPadDownSimHits_v12_trackIDContribs_,
                           &b_TriggerPadDownSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.incidentIDContribs_",
                           TriggerPadDownSimHits_v12_incidentIDContribs_,
                           &b_TriggerPadDownSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.pdgCodeContribs_",
                           TriggerPadDownSimHits_v12_pdgCodeContribs_,
                           &b_TriggerPadDownSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.edepContribs_",
                           TriggerPadDownSimHits_v12_edepContribs_,
                           &b_TriggerPadDownSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.timeContribs_",
                           TriggerPadDownSimHits_v12_timeContribs_,
                           &b_TriggerPadDownSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("TriggerPadDownSimHits_v12.nContribs_",
                           TriggerPadDownSimHits_v12_nContribs_,
                           &b_TriggerPadDownSimHits_v12_nContribs_);
  fChain->SetBranchAddress("RecoilSimHits_v12", &RecoilSimHits_v12_,
                           &b_RecoilSimHits_v12_);
  fChain->SetBranchAddress("RecoilSimHits_v12.id_", RecoilSimHits_v12_id_,
                           &b_RecoilSimHits_v12_id_);
  fChain->SetBranchAddress("RecoilSimHits_v12.layerID_",
                           RecoilSimHits_v12_layerID_,
                           &b_RecoilSimHits_v12_layerID_);
  fChain->SetBranchAddress("RecoilSimHits_v12.moduleID_",
                           RecoilSimHits_v12_moduleID_,
                           &b_RecoilSimHits_v12_moduleID_);
  fChain->SetBranchAddress("RecoilSimHits_v12.edep_", RecoilSimHits_v12_edep_,
                           &b_RecoilSimHits_v12_edep_);
  fChain->SetBranchAddress("RecoilSimHits_v12.time_", RecoilSimHits_v12_time_,
                           &b_RecoilSimHits_v12_time_);
  fChain->SetBranchAddress("RecoilSimHits_v12.px_", RecoilSimHits_v12_px_,
                           &b_RecoilSimHits_v12_px_);
  fChain->SetBranchAddress("RecoilSimHits_v12.py_", RecoilSimHits_v12_py_,
                           &b_RecoilSimHits_v12_py_);
  fChain->SetBranchAddress("RecoilSimHits_v12.pz_", RecoilSimHits_v12_pz_,
                           &b_RecoilSimHits_v12_pz_);
  fChain->SetBranchAddress("RecoilSimHits_v12.energy_",
                           RecoilSimHits_v12_energy_,
                           &b_RecoilSimHits_v12_energy_);
  fChain->SetBranchAddress("RecoilSimHits_v12.x_", RecoilSimHits_v12_x_,
                           &b_RecoilSimHits_v12_x_);
  fChain->SetBranchAddress("RecoilSimHits_v12.y_", RecoilSimHits_v12_y_,
                           &b_RecoilSimHits_v12_y_);
  fChain->SetBranchAddress("RecoilSimHits_v12.z_", RecoilSimHits_v12_z_,
                           &b_RecoilSimHits_v12_z_);
  fChain->SetBranchAddress("RecoilSimHits_v12.pathLength_",
                           RecoilSimHits_v12_pathLength_,
                           &b_RecoilSimHits_v12_pathLength_);
  fChain->SetBranchAddress(
      "RecoilSimHits_v12.trackID_", RecoilSimHits_v12_trackID_,
      &b_RecoilSimHits_v12_trackID_);  // the info of what particle caused this
                                       // hit??
  fChain->SetBranchAddress("RecoilSimHits_v12.pdgID_", RecoilSimHits_v12_pdgID_,
                           &b_RecoilSimHits_v12_pdgID_);
  fChain->SetBranchAddress("EcalSimHits_v12", &EcalSimHits_v12_,
                           &b_EcalSimHits_v12_);
  fChain->SetBranchAddress("EcalSimHits_v12.id_", EcalSimHits_v12_id_,
                           &b_EcalSimHits_v12_id_);
  fChain->SetBranchAddress("EcalSimHits_v12.edep_", EcalSimHits_v12_edep_,
                           &b_EcalSimHits_v12_edep_);
  fChain->SetBranchAddress("EcalSimHits_v12.x_", EcalSimHits_v12_x_,
                           &b_EcalSimHits_v12_x_);
  fChain->SetBranchAddress("EcalSimHits_v12.y_", EcalSimHits_v12_y_,
                           &b_EcalSimHits_v12_y_);
  fChain->SetBranchAddress("EcalSimHits_v12.z_", EcalSimHits_v12_z_,
                           &b_EcalSimHits_v12_z_);
  fChain->SetBranchAddress("EcalSimHits_v12.time_", EcalSimHits_v12_time_,
                           &b_EcalSimHits_v12_time_);
  fChain->SetBranchAddress("EcalSimHits_v12.trackIDContribs_",
                           EcalSimHits_v12_trackIDContribs_,
                           &b_EcalSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("EcalSimHits_v12.incidentIDContribs_",
                           EcalSimHits_v12_incidentIDContribs_,
                           &b_EcalSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("EcalSimHits_v12.pdgCodeContribs_",
                           EcalSimHits_v12_pdgCodeContribs_,
                           &b_EcalSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("EcalSimHits_v12.edepContribs_",
                           EcalSimHits_v12_edepContribs_,
                           &b_EcalSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("EcalSimHits_v12.timeContribs_",
                           EcalSimHits_v12_timeContribs_,
                           &b_EcalSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("EcalSimHits_v12.nContribs_",
                           EcalSimHits_v12_nContribs_,
                           &b_EcalSimHits_v12_nContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12", &HcalSimHits_v12_,
                           &b_HcalSimHits_v12_);
  fChain->SetBranchAddress("HcalSimHits_v12.id_", HcalSimHits_v12_id_,
                           &b_HcalSimHits_v12_id_);
  fChain->SetBranchAddress("HcalSimHits_v12.edep_", HcalSimHits_v12_edep_,
                           &b_HcalSimHits_v12_edep_);
  fChain->SetBranchAddress("HcalSimHits_v12.x_", HcalSimHits_v12_x_,
                           &b_HcalSimHits_v12_x_);
  fChain->SetBranchAddress("HcalSimHits_v12.y_", HcalSimHits_v12_y_,
                           &b_HcalSimHits_v12_y_);
  fChain->SetBranchAddress("HcalSimHits_v12.z_", HcalSimHits_v12_z_,
                           &b_HcalSimHits_v12_z_);
  fChain->SetBranchAddress("HcalSimHits_v12.time_", HcalSimHits_v12_time_,
                           &b_HcalSimHits_v12_time_);
  fChain->SetBranchAddress("HcalSimHits_v12.trackIDContribs_",
                           HcalSimHits_v12_trackIDContribs_,
                           &b_HcalSimHits_v12_trackIDContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12.incidentIDContribs_",
                           HcalSimHits_v12_incidentIDContribs_,
                           &b_HcalSimHits_v12_incidentIDContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12.pdgCodeContribs_",
                           HcalSimHits_v12_pdgCodeContribs_,
                           &b_HcalSimHits_v12_pdgCodeContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12.edepContribs_",
                           HcalSimHits_v12_edepContribs_,
                           &b_HcalSimHits_v12_edepContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12.timeContribs_",
                           HcalSimHits_v12_timeContribs_,
                           &b_HcalSimHits_v12_timeContribs_);
  fChain->SetBranchAddress("HcalSimHits_v12.nContribs_",
                           HcalSimHits_v12_nContribs_,
                           &b_HcalSimHits_v12_nContribs_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12",
                           &MagnetScoringPlaneHits_v12_,
                           &b_MagnetScoringPlaneHits_v12_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.id_",
                           MagnetScoringPlaneHits_v12_id_,
                           &b_MagnetScoringPlaneHits_v12_id_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.layerID_",
                           MagnetScoringPlaneHits_v12_layerID_,
                           &b_MagnetScoringPlaneHits_v12_layerID_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.moduleID_",
                           MagnetScoringPlaneHits_v12_moduleID_,
                           &b_MagnetScoringPlaneHits_v12_moduleID_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.edep_",
                           MagnetScoringPlaneHits_v12_edep_,
                           &b_MagnetScoringPlaneHits_v12_edep_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.time_",
                           MagnetScoringPlaneHits_v12_time_,
                           &b_MagnetScoringPlaneHits_v12_time_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.px_",
                           MagnetScoringPlaneHits_v12_px_,
                           &b_MagnetScoringPlaneHits_v12_px_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.py_",
                           MagnetScoringPlaneHits_v12_py_,
                           &b_MagnetScoringPlaneHits_v12_py_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.pz_",
                           MagnetScoringPlaneHits_v12_pz_,
                           &b_MagnetScoringPlaneHits_v12_pz_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.energy_",
                           MagnetScoringPlaneHits_v12_energy_,
                           &b_MagnetScoringPlaneHits_v12_energy_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.x_",
                           MagnetScoringPlaneHits_v12_x_,
                           &b_MagnetScoringPlaneHits_v12_x_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.y_",
                           MagnetScoringPlaneHits_v12_y_,
                           &b_MagnetScoringPlaneHits_v12_y_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.z_",
                           MagnetScoringPlaneHits_v12_z_,
                           &b_MagnetScoringPlaneHits_v12_z_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.pathLength_",
                           MagnetScoringPlaneHits_v12_pathLength_,
                           &b_MagnetScoringPlaneHits_v12_pathLength_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.trackID_",
                           MagnetScoringPlaneHits_v12_trackID_,
                           &b_MagnetScoringPlaneHits_v12_trackID_);
  fChain->SetBranchAddress("MagnetScoringPlaneHits_v12.pdgID_",
                           MagnetScoringPlaneHits_v12_pdgID_,
                           &b_MagnetScoringPlaneHits_v12_pdgID_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12",
                           &TargetScoringPlaneHits_v12_,
                           &b_TargetScoringPlaneHits_v12_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.id_",
                           TargetScoringPlaneHits_v12_id_,
                           &b_TargetScoringPlaneHits_v12_id_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.layerID_",
                           TargetScoringPlaneHits_v12_layerID_,
                           &b_TargetScoringPlaneHits_v12_layerID_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.moduleID_",
                           TargetScoringPlaneHits_v12_moduleID_,
                           &b_TargetScoringPlaneHits_v12_moduleID_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.edep_",
                           TargetScoringPlaneHits_v12_edep_,
                           &b_TargetScoringPlaneHits_v12_edep_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.time_",
                           TargetScoringPlaneHits_v12_time_,
                           &b_TargetScoringPlaneHits_v12_time_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.px_",
                           TargetScoringPlaneHits_v12_px_,
                           &b_TargetScoringPlaneHits_v12_px_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.py_",
                           TargetScoringPlaneHits_v12_py_,
                           &b_TargetScoringPlaneHits_v12_py_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.pz_",
                           TargetScoringPlaneHits_v12_pz_,
                           &b_TargetScoringPlaneHits_v12_pz_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.energy_",
                           TargetScoringPlaneHits_v12_energy_,
                           &b_TargetScoringPlaneHits_v12_energy_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.x_",
                           TargetScoringPlaneHits_v12_x_,
                           &b_TargetScoringPlaneHits_v12_x_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.y_",
                           TargetScoringPlaneHits_v12_y_,
                           &b_TargetScoringPlaneHits_v12_y_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.z_",
                           TargetScoringPlaneHits_v12_z_,
                           &b_TargetScoringPlaneHits_v12_z_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.pathLength_",
                           TargetScoringPlaneHits_v12_pathLength_,
                           &b_TargetScoringPlaneHits_v12_pathLength_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.trackID_",
                           TargetScoringPlaneHits_v12_trackID_,
                           &b_TargetScoringPlaneHits_v12_trackID_);
  fChain->SetBranchAddress("TargetScoringPlaneHits_v12.pdgID_",
                           TargetScoringPlaneHits_v12_pdgID_,
                           &b_TargetScoringPlaneHits_v12_pdgID_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12",
                           &TrackerScoringPlaneHits_v12_,
                           &b_TrackerScoringPlaneHits_v12_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.id_",
                           TrackerScoringPlaneHits_v12_id_,
                           &b_TrackerScoringPlaneHits_v12_id_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.layerID_",
                           TrackerScoringPlaneHits_v12_layerID_,
                           &b_TrackerScoringPlaneHits_v12_layerID_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.moduleID_",
                           TrackerScoringPlaneHits_v12_moduleID_,
                           &b_TrackerScoringPlaneHits_v12_moduleID_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.edep_",
                           TrackerScoringPlaneHits_v12_edep_,
                           &b_TrackerScoringPlaneHits_v12_edep_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.time_",
                           TrackerScoringPlaneHits_v12_time_,
                           &b_TrackerScoringPlaneHits_v12_time_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.px_",
                           TrackerScoringPlaneHits_v12_px_,
                           &b_TrackerScoringPlaneHits_v12_px_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.py_",
                           TrackerScoringPlaneHits_v12_py_,
                           &b_TrackerScoringPlaneHits_v12_py_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.pz_",
                           TrackerScoringPlaneHits_v12_pz_,
                           &b_TrackerScoringPlaneHits_v12_pz_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.energy_",
                           TrackerScoringPlaneHits_v12_energy_,
                           &b_TrackerScoringPlaneHits_v12_energy_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.x_",
                           TrackerScoringPlaneHits_v12_x_,
                           &b_TrackerScoringPlaneHits_v12_x_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.y_",
                           TrackerScoringPlaneHits_v12_y_,
                           &b_TrackerScoringPlaneHits_v12_y_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.z_",
                           TrackerScoringPlaneHits_v12_z_,
                           &b_TrackerScoringPlaneHits_v12_z_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.pathLength_",
                           TrackerScoringPlaneHits_v12_pathLength_,
                           &b_TrackerScoringPlaneHits_v12_pathLength_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.trackID_",
                           TrackerScoringPlaneHits_v12_trackID_,
                           &b_TrackerScoringPlaneHits_v12_trackID_);
  fChain->SetBranchAddress("TrackerScoringPlaneHits_v12.pdgID_",
                           TrackerScoringPlaneHits_v12_pdgID_,
                           &b_TrackerScoringPlaneHits_v12_pdgID_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12",
                           &EcalScoringPlaneHits_v12_,
                           &b_EcalScoringPlaneHits_v12_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.id_",
                           EcalScoringPlaneHits_v12_id_,
                           &b_EcalScoringPlaneHits_v12_id_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.layerID_",
                           EcalScoringPlaneHits_v12_layerID_,
                           &b_EcalScoringPlaneHits_v12_layerID_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.moduleID_",
                           EcalScoringPlaneHits_v12_moduleID_,
                           &b_EcalScoringPlaneHits_v12_moduleID_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.edep_",
                           EcalScoringPlaneHits_v12_edep_,
                           &b_EcalScoringPlaneHits_v12_edep_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.time_",
                           EcalScoringPlaneHits_v12_time_,
                           &b_EcalScoringPlaneHits_v12_time_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.px_",
                           EcalScoringPlaneHits_v12_px_,
                           &b_EcalScoringPlaneHits_v12_px_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.py_",
                           EcalScoringPlaneHits_v12_py_,
                           &b_EcalScoringPlaneHits_v12_py_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.pz_",
                           EcalScoringPlaneHits_v12_pz_,
                           &b_EcalScoringPlaneHits_v12_pz_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.energy_",
                           EcalScoringPlaneHits_v12_energy_,
                           &b_EcalScoringPlaneHits_v12_energy_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.x_",
                           EcalScoringPlaneHits_v12_x_,
                           &b_EcalScoringPlaneHits_v12_x_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.y_",
                           EcalScoringPlaneHits_v12_y_,
                           &b_EcalScoringPlaneHits_v12_y_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.z_",
                           EcalScoringPlaneHits_v12_z_,
                           &b_EcalScoringPlaneHits_v12_z_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.pathLength_",
                           EcalScoringPlaneHits_v12_pathLength_,
                           &b_EcalScoringPlaneHits_v12_pathLength_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.trackID_",
                           EcalScoringPlaneHits_v12_trackID_,
                           &b_EcalScoringPlaneHits_v12_trackID_);
  fChain->SetBranchAddress("EcalScoringPlaneHits_v12.pdgID_",
                           EcalScoringPlaneHits_v12_pdgID_,
                           &b_EcalScoringPlaneHits_v12_pdgID_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12",
                           &HcalScoringPlaneHits_v12_,
                           &b_HcalScoringPlaneHits_v12_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.id_",
                           HcalScoringPlaneHits_v12_id_,
                           &b_HcalScoringPlaneHits_v12_id_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.layerID_",
                           HcalScoringPlaneHits_v12_layerID_,
                           &b_HcalScoringPlaneHits_v12_layerID_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.moduleID_",
                           HcalScoringPlaneHits_v12_moduleID_,
                           &b_HcalScoringPlaneHits_v12_moduleID_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.edep_",
                           HcalScoringPlaneHits_v12_edep_,
                           &b_HcalScoringPlaneHits_v12_edep_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.time_",
                           HcalScoringPlaneHits_v12_time_,
                           &b_HcalScoringPlaneHits_v12_time_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.px_",
                           HcalScoringPlaneHits_v12_px_,
                           &b_HcalScoringPlaneHits_v12_px_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.py_",
                           HcalScoringPlaneHits_v12_py_,
                           &b_HcalScoringPlaneHits_v12_py_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.pz_",
                           HcalScoringPlaneHits_v12_pz_,
                           &b_HcalScoringPlaneHits_v12_pz_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.energy_",
                           HcalScoringPlaneHits_v12_energy_,
                           &b_HcalScoringPlaneHits_v12_energy_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.x_",
                           HcalScoringPlaneHits_v12_x_,
                           &b_HcalScoringPlaneHits_v12_x_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.y_",
                           HcalScoringPlaneHits_v12_y_,
                           &b_HcalScoringPlaneHits_v12_y_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.z_",
                           HcalScoringPlaneHits_v12_z_,
                           &b_HcalScoringPlaneHits_v12_z_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.pathLength_",
                           HcalScoringPlaneHits_v12_pathLength_,
                           &b_HcalScoringPlaneHits_v12_pathLength_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.trackID_",
                           HcalScoringPlaneHits_v12_trackID_,
                           &b_HcalScoringPlaneHits_v12_trackID_);
  fChain->SetBranchAddress("HcalScoringPlaneHits_v12.pdgID_",
                           HcalScoringPlaneHits_v12_pdgID_,
                           &b_HcalScoringPlaneHits_v12_pdgID_);
  fChain->SetBranchAddress("eventNumber_", &eventNumber_,
                           &b_EventHeader_eventNumber_);
  fChain->SetBranchAddress("run_", &run_, &b_EventHeader_run_);
  fChain->SetBranchAddress("timestamp_.fSec", &timestamp__fSec,
                           &b_EventHeader_timestamp__fSec);
  fChain->SetBranchAddress("timestamp_.fNanoSec", &timestamp__fNanoSec,
                           &b_EventHeader_timestamp__fNanoSec);
  fChain->SetBranchAddress("weight_", &weight_, &b_EventHeader_weight_);
  fChain->SetBranchAddress("isRealData_", &isRealData_,
                           &b_EventHeader_isRealData_);
  fChain->SetBranchAddress("intParameters_", &intParameters__,
                           &b_EventHeader_intParameters__);
  fChain->SetBranchAddress("intParameters_.first", &intParameters__first,
                           &b_intParameters__first);
  fChain->SetBranchAddress("intParameters_.second", &intParameters__second,
                           &b_intParameters__second);
  fChain->SetBranchAddress("floatParameters_", &floatParameters__,
                           &b_EventHeader_floatParameters__);
  fChain->SetBranchAddress("floatParameters_.first", &floatParameters__first,
                           &b_floatParameters__first);
  fChain->SetBranchAddress("floatParameters_.second", &floatParameters__second,
                           &b_floatParameters__second);
  fChain->SetBranchAddress("stringParameters_", &stringParameters__,
                           &b_EventHeader_stringParameters__);
  fChain->SetBranchAddress("stringParameters_.first", stringParameters__first,
                           &b_stringParameters__first);
  fChain->SetBranchAddress("stringParameters_.second", stringParameters__second,
                           &b_stringParameters__second);
  fChain->SetBranchAddress("channelIDs_", &channelIDs_,
                           &b_EcalDigis_reco_channelIDs_);
  fChain->SetBranchAddress("samples_", &samples__, &b_EcalDigis_reco_samples__);
  fChain->SetBranchAddress("samples_.word_", samples__word_, &b_samples__word_);
  fChain->SetBranchAddress("numSamplesPerDigi_", &numSamplesPerDigi_,
                           &b_EcalDigis_reco_numSamplesPerDigi_);
  fChain->SetBranchAddress("sampleOfInterest_", &sampleOfInterest_,
                           &b_EcalDigis_reco_sampleOfInterest_);
  fChain->SetBranchAddress("EcalRecHits_reco", &EcalRecHits_reco_,
                           &b_EcalRecHits_reco_);
  fChain->SetBranchAddress("EcalRecHits_reco.id_", EcalRecHits_reco_id_,
                           &b_EcalRecHits_reco_id_);
  fChain->SetBranchAddress("EcalRecHits_reco.amplitude_",
                           EcalRecHits_reco_amplitude_,
                           &b_EcalRecHits_reco_amplitude_);
  fChain->SetBranchAddress("EcalRecHits_reco.energy_", EcalRecHits_reco_energy_,
                           &b_EcalRecHits_reco_energy_);
  fChain->SetBranchAddress("EcalRecHits_reco.time_", EcalRecHits_reco_time_,
                           &b_EcalRecHits_reco_time_);
  fChain->SetBranchAddress("EcalRecHits_reco.xpos_", EcalRecHits_reco_xpos_,
                           &b_EcalRecHits_reco_xpos_);
  fChain->SetBranchAddress("EcalRecHits_reco.ypos_", EcalRecHits_reco_ypos_,
                           &b_EcalRecHits_reco_ypos_);
  fChain->SetBranchAddress("EcalRecHits_reco.zpos_", EcalRecHits_reco_zpos_,
                           &b_EcalRecHits_reco_zpos_);
  fChain->SetBranchAddress("EcalRecHits_reco.isNoise_",
                           EcalRecHits_reco_isNoise_,
                           &b_EcalRecHits_reco_isNoise_);
  fChain->SetBranchAddress("passesVeto_", &passesVeto_,
                           &b_EcalVeto_reco_passesVeto_);
  fChain->SetBranchAddress("nReadoutHits_", &nReadoutHits_,
                           &b_EcalVeto_reco_nReadoutHits_);
  fChain->SetBranchAddress("deepestLayerHit_", &deepestLayerHit_,
                           &b_EcalVeto_reco_deepestLayerHit_);
  fChain->SetBranchAddress("summedDet_", &summedDet_,
                           &b_EcalVeto_reco_summedDet_);
  fChain->SetBranchAddress("summedTightIso_", &summedTightIso_,
                           &b_EcalVeto_reco_summedTightIso_);
  fChain->SetBranchAddress("maxCellDep_", &maxCellDep_,
                           &b_EcalVeto_reco_maxCellDep_);
  fChain->SetBranchAddress("showerRMS_", &showerRMS_,
                           &b_EcalVeto_reco_showerRMS_);
  fChain->SetBranchAddress("xStd_", &xStd_, &b_EcalVeto_reco_xStd_);
  fChain->SetBranchAddress("yStd_", &yStd_, &b_EcalVeto_reco_yStd_);
  fChain->SetBranchAddress("avgLayerHit_", &avgLayerHit_,
                           &b_EcalVeto_reco_avgLayerHit_);
  fChain->SetBranchAddress("stdLayerHit_", &stdLayerHit_,
                           &b_EcalVeto_reco_stdLayerHit_);
  fChain->SetBranchAddress("ecalBackEnergy_", &ecalBackEnergy_,
                           &b_EcalVeto_reco_ecalBackEnergy_);
  fChain->SetBranchAddress("electronContainmentEnergy_",
                           &electronContainmentEnergy_,
                           &b_EcalVeto_reco_electronContainmentEnergy_);
  fChain->SetBranchAddress("photonContainmentEnergy_",
                           &photonContainmentEnergy_,
                           &b_EcalVeto_reco_photonContainmentEnergy_);
  fChain->SetBranchAddress("outsideContainmentEnergy_",
                           &outsideContainmentEnergy_,
                           &b_EcalVeto_reco_outsideContainmentEnergy_);
  fChain->SetBranchAddress("outsideContainmentNHits_",
                           &outsideContainmentNHits_,
                           &b_EcalVeto_reco_outsideContainmentNHits_);
  fChain->SetBranchAddress("outsideContainmentXStd_", &outsideContainmentXStd_,
                           &b_EcalVeto_reco_outsideContainmentXStd_);
  fChain->SetBranchAddress("outsideContainmentYStd_", &outsideContainmentYStd_,
                           &b_EcalVeto_reco_outsideContainmentYStd_);
  fChain->SetBranchAddress("discValue_", &discValue_,
                           &b_EcalVeto_reco_discValue_);
  fChain->SetBranchAddress("recoilPx_", &recoilPx_, &b_EcalVeto_reco_recoilPx_);
  fChain->SetBranchAddress("recoilPy_", &recoilPy_, &b_EcalVeto_reco_recoilPy_);
  fChain->SetBranchAddress("recoilPz_", &recoilPz_, &b_EcalVeto_reco_recoilPz_);
  fChain->SetBranchAddress("recoilX_", &recoilX_, &b_EcalVeto_reco_recoilX_);
  fChain->SetBranchAddress("recoilY_", &recoilY_, &b_EcalVeto_reco_recoilY_);
  fChain->SetBranchAddress("ecalLayerEdepReadout_", &ecalLayerEdepReadout_,
                           &b_EcalVeto_reco_ecalLayerEdepReadout_);
  fChain->SetBranchAddress("HcalRecHits_reco", &HcalRecHits_reco_,
                           &b_HcalRecHits_reco_);
  fChain->SetBranchAddress("HcalRecHits_reco.id_", HcalRecHits_reco_id_,
                           &b_HcalRecHits_reco_id_);
  fChain->SetBranchAddress("HcalRecHits_reco.amplitude_",
                           HcalRecHits_reco_amplitude_,
                           &b_HcalRecHits_reco_amplitude_);
  fChain->SetBranchAddress("HcalRecHits_reco.energy_", HcalRecHits_reco_energy_,
                           &b_HcalRecHits_reco_energy_);
  fChain->SetBranchAddress("HcalRecHits_reco.time_", HcalRecHits_reco_time_,
                           &b_HcalRecHits_reco_time_);
  fChain->SetBranchAddress("HcalRecHits_reco.xpos_", HcalRecHits_reco_xpos_,
                           &b_HcalRecHits_reco_xpos_);
  fChain->SetBranchAddress("HcalRecHits_reco.ypos_", HcalRecHits_reco_ypos_,
                           &b_HcalRecHits_reco_ypos_);
  fChain->SetBranchAddress("HcalRecHits_reco.zpos_", HcalRecHits_reco_zpos_,
                           &b_HcalRecHits_reco_zpos_);
  fChain->SetBranchAddress("HcalRecHits_reco.isNoise_",
                           HcalRecHits_reco_isNoise_,
                           &b_HcalRecHits_reco_isNoise_);
  fChain->SetBranchAddress("HcalRecHits_reco.pe_", HcalRecHits_reco_pe_,
                           &b_HcalRecHits_reco_pe_);
  fChain->SetBranchAddress("HcalRecHits_reco.minpe_", HcalRecHits_reco_minpe_,
                           &b_HcalRecHits_reco_minpe_);
  fChain->SetBranchAddress("maxPEHit_.id_", &maxPEHit__id_,
                           &b_HcalVeto_reco_maxPEHit__id_);
  fChain->SetBranchAddress("maxPEHit_.amplitude_", &maxPEHit__amplitude_,
                           &b_HcalVeto_reco_maxPEHit__amplitude_);
  fChain->SetBranchAddress("maxPEHit_.energy_", &maxPEHit__energy_,
                           &b_HcalVeto_reco_maxPEHit__energy_);
  fChain->SetBranchAddress("maxPEHit_.time_", &maxPEHit__time_,
                           &b_HcalVeto_reco_maxPEHit__time_);
  fChain->SetBranchAddress("maxPEHit_.xpos_", &maxPEHit__xpos_,
                           &b_HcalVeto_reco_maxPEHit__xpos_);
  fChain->SetBranchAddress("maxPEHit_.ypos_", &maxPEHit__ypos_,
                           &b_HcalVeto_reco_maxPEHit__ypos_);
  fChain->SetBranchAddress("maxPEHit_.zpos_", &maxPEHit__zpos_,
                           &b_HcalVeto_reco_maxPEHit__zpos_);
  fChain->SetBranchAddress("maxPEHit_.isNoise_", &maxPEHit__isNoise_,
                           &b_HcalVeto_reco_maxPEHit__isNoise_);
  fChain->SetBranchAddress("maxPEHit_.pe_", &maxPEHit__pe_,
                           &b_HcalVeto_reco_maxPEHit__pe_);
  fChain->SetBranchAddress("maxPEHit_.minpe_", &maxPEHit__minpe_,
                           &b_HcalVeto_reco_maxPEHit__minpe_);
  //    fChain->SetBranchAddress("passesVeto_", &passesVeto_,
  //    &b_HcalVeto_reco_passesVeto_);
  fChain->SetBranchAddress("trigScintDigisTag_reco", &trigScintDigisTag_reco_,
                           &b_trigScintDigisTag_reco_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.id_",
                           trigScintDigisTag_reco_id_,
                           &b_trigScintDigisTag_reco_id_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.amplitude_",
                           trigScintDigisTag_reco_amplitude_,
                           &b_trigScintDigisTag_reco_amplitude_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.energy_",
                           trigScintDigisTag_reco_energy_,
                           &b_trigScintDigisTag_reco_energy_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.time_",
                           trigScintDigisTag_reco_time_,
                           &b_trigScintDigisTag_reco_time_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.xpos_",
                           trigScintDigisTag_reco_xpos_,
                           &b_trigScintDigisTag_reco_xpos_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.ypos_",
                           trigScintDigisTag_reco_ypos_,
                           &b_trigScintDigisTag_reco_ypos_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.zpos_",
                           trigScintDigisTag_reco_zpos_,
                           &b_trigScintDigisTag_reco_zpos_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.isNoise_",
                           trigScintDigisTag_reco_isNoise_,
                           &b_trigScintDigisTag_reco_isNoise_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.pe_",
                           trigScintDigisTag_reco_pe_,
                           &b_trigScintDigisTag_reco_pe_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.minpe_",
                           trigScintDigisTag_reco_minpe_,
                           &b_trigScintDigisTag_reco_minpe_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.barID_",
                           trigScintDigisTag_reco_barID_,
                           &b_trigScintDigisTag_reco_barID_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.moduleID_",
                           trigScintDigisTag_reco_moduleID_,
                           &b_trigScintDigisTag_reco_moduleID_);
  fChain->SetBranchAddress("trigScintDigisTag_reco.beamEfrac_",
                           trigScintDigisTag_reco_beamEfrac_,
                           &b_trigScintDigisTag_reco_beamEfrac_);
  fChain->SetBranchAddress("trigScintDigisUp_reco", &trigScintDigisUp_reco_,
                           &b_trigScintDigisUp_reco_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.id_",
                           trigScintDigisUp_reco_id_,
                           &b_trigScintDigisUp_reco_id_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.amplitude_",
                           trigScintDigisUp_reco_amplitude_,
                           &b_trigScintDigisUp_reco_amplitude_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.energy_",
                           trigScintDigisUp_reco_energy_,
                           &b_trigScintDigisUp_reco_energy_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.time_",
                           trigScintDigisUp_reco_time_,
                           &b_trigScintDigisUp_reco_time_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.xpos_",
                           trigScintDigisUp_reco_xpos_,
                           &b_trigScintDigisUp_reco_xpos_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.ypos_",
                           trigScintDigisUp_reco_ypos_,
                           &b_trigScintDigisUp_reco_ypos_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.zpos_",
                           trigScintDigisUp_reco_zpos_,
                           &b_trigScintDigisUp_reco_zpos_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.isNoise_",
                           trigScintDigisUp_reco_isNoise_,
                           &b_trigScintDigisUp_reco_isNoise_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.pe_",
                           trigScintDigisUp_reco_pe_,
                           &b_trigScintDigisUp_reco_pe_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.minpe_",
                           trigScintDigisUp_reco_minpe_,
                           &b_trigScintDigisUp_reco_minpe_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.barID_",
                           trigScintDigisUp_reco_barID_,
                           &b_trigScintDigisUp_reco_barID_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.moduleID_",
                           trigScintDigisUp_reco_moduleID_,
                           &b_trigScintDigisUp_reco_moduleID_);
  fChain->SetBranchAddress("trigScintDigisUp_reco.beamEfrac_",
                           trigScintDigisUp_reco_beamEfrac_,
                           &b_trigScintDigisUp_reco_beamEfrac_);
  fChain->SetBranchAddress("trigScintDigisDn_reco", &trigScintDigisDn_reco_,
                           &b_trigScintDigisDn_reco_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.id_",
                           trigScintDigisDn_reco_id_,
                           &b_trigScintDigisDn_reco_id_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.amplitude_",
                           trigScintDigisDn_reco_amplitude_,
                           &b_trigScintDigisDn_reco_amplitude_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.energy_",
                           trigScintDigisDn_reco_energy_,
                           &b_trigScintDigisDn_reco_energy_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.time_",
                           trigScintDigisDn_reco_time_,
                           &b_trigScintDigisDn_reco_time_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.xpos_",
                           trigScintDigisDn_reco_xpos_,
                           &b_trigScintDigisDn_reco_xpos_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.ypos_",
                           trigScintDigisDn_reco_ypos_,
                           &b_trigScintDigisDn_reco_ypos_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.zpos_",
                           trigScintDigisDn_reco_zpos_,
                           &b_trigScintDigisDn_reco_zpos_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.isNoise_",
                           trigScintDigisDn_reco_isNoise_,
                           &b_trigScintDigisDn_reco_isNoise_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.pe_",
                           trigScintDigisDn_reco_pe_,
                           &b_trigScintDigisDn_reco_pe_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.minpe_",
                           trigScintDigisDn_reco_minpe_,
                           &b_trigScintDigisDn_reco_minpe_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.barID_",
                           trigScintDigisDn_reco_barID_,
                           &b_trigScintDigisDn_reco_barID_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.moduleID_",
                           trigScintDigisDn_reco_moduleID_,
                           &b_trigScintDigisDn_reco_moduleID_);
  fChain->SetBranchAddress("trigScintDigisDn_reco.beamEfrac_",
                           trigScintDigisDn_reco_beamEfrac_,
                           &b_trigScintDigisDn_reco_beamEfrac_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco",
                           &TriggerPadTaggerClusters_reco_,
                           &b_TriggerPadTaggerClusters_reco_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.hitIDs_",
                           TriggerPadTaggerClusters_reco_hitIDs_,
                           &b_TriggerPadTaggerClusters_reco_hitIDs_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.energy_",
                           TriggerPadTaggerClusters_reco_energy_,
                           &b_TriggerPadTaggerClusters_reco_energy_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.nHits_",
                           TriggerPadTaggerClusters_reco_nHits_,
                           &b_TriggerPadTaggerClusters_reco_nHits_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.PE_",
                           TriggerPadTaggerClusters_reco_PE_,
                           &b_TriggerPadTaggerClusters_reco_PE_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.seed_",
                           TriggerPadTaggerClusters_reco_seed_,
                           &b_TriggerPadTaggerClusters_reco_seed_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.centroid_",
                           TriggerPadTaggerClusters_reco_centroid_,
                           &b_TriggerPadTaggerClusters_reco_centroid_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.centroidX_",
                           TriggerPadTaggerClusters_reco_centroidX_,
                           &b_TriggerPadTaggerClusters_reco_centroidX_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.centroidY_",
                           TriggerPadTaggerClusters_reco_centroidY_,
                           &b_TriggerPadTaggerClusters_reco_centroidY_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.centroidZ_",
                           TriggerPadTaggerClusters_reco_centroidZ_,
                           &b_TriggerPadTaggerClusters_reco_centroidZ_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.beamEfrac_",
                           TriggerPadTaggerClusters_reco_beamEfrac_,
                           &b_TriggerPadTaggerClusters_reco_beamEfrac_);
  fChain->SetBranchAddress("TriggerPadTaggerClusters_reco.time_",
                           TriggerPadTaggerClusters_reco_time_,
                           &b_TriggerPadTaggerClusters_reco_time_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco",
                           &TriggerPadUpClusters_reco_,
                           &b_TriggerPadUpClusters_reco_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.hitIDs_",
                           TriggerPadUpClusters_reco_hitIDs_,
                           &b_TriggerPadUpClusters_reco_hitIDs_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.energy_",
                           TriggerPadUpClusters_reco_energy_,
                           &b_TriggerPadUpClusters_reco_energy_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.nHits_",
                           TriggerPadUpClusters_reco_nHits_,
                           &b_TriggerPadUpClusters_reco_nHits_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.PE_",
                           TriggerPadUpClusters_reco_PE_,
                           &b_TriggerPadUpClusters_reco_PE_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.seed_",
                           TriggerPadUpClusters_reco_seed_,
                           &b_TriggerPadUpClusters_reco_seed_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.centroid_",
                           TriggerPadUpClusters_reco_centroid_,
                           &b_TriggerPadUpClusters_reco_centroid_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.centroidX_",
                           TriggerPadUpClusters_reco_centroidX_,
                           &b_TriggerPadUpClusters_reco_centroidX_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.centroidY_",
                           TriggerPadUpClusters_reco_centroidY_,
                           &b_TriggerPadUpClusters_reco_centroidY_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.centroidZ_",
                           TriggerPadUpClusters_reco_centroidZ_,
                           &b_TriggerPadUpClusters_reco_centroidZ_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.beamEfrac_",
                           TriggerPadUpClusters_reco_beamEfrac_,
                           &b_TriggerPadUpClusters_reco_beamEfrac_);
  fChain->SetBranchAddress("TriggerPadUpClusters_reco.time_",
                           TriggerPadUpClusters_reco_time_,
                           &b_TriggerPadUpClusters_reco_time_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco",
                           &TriggerPadDownClusters_reco_,
                           &b_TriggerPadDownClusters_reco_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.hitIDs_",
                           TriggerPadDownClusters_reco_hitIDs_,
                           &b_TriggerPadDownClusters_reco_hitIDs_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.energy_",
                           TriggerPadDownClusters_reco_energy_,
                           &b_TriggerPadDownClusters_reco_energy_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.nHits_",
                           TriggerPadDownClusters_reco_nHits_,
                           &b_TriggerPadDownClusters_reco_nHits_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.PE_",
                           TriggerPadDownClusters_reco_PE_,
                           &b_TriggerPadDownClusters_reco_PE_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.seed_",
                           TriggerPadDownClusters_reco_seed_,
                           &b_TriggerPadDownClusters_reco_seed_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.centroid_",
                           TriggerPadDownClusters_reco_centroid_,
                           &b_TriggerPadDownClusters_reco_centroid_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.centroidX_",
                           TriggerPadDownClusters_reco_centroidX_,
                           &b_TriggerPadDownClusters_reco_centroidX_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.centroidY_",
                           TriggerPadDownClusters_reco_centroidY_,
                           &b_TriggerPadDownClusters_reco_centroidY_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.centroidZ_",
                           TriggerPadDownClusters_reco_centroidZ_,
                           &b_TriggerPadDownClusters_reco_centroidZ_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.beamEfrac_",
                           TriggerPadDownClusters_reco_beamEfrac_,
                           &b_TriggerPadDownClusters_reco_beamEfrac_);
  fChain->SetBranchAddress("TriggerPadDownClusters_reco.time_",
                           TriggerPadDownClusters_reco_time_,
                           &b_TriggerPadDownClusters_reco_time_);
  fChain->SetBranchAddress("TriggerPadTracks_reco", &TriggerPadTracks_reco_,
                           &b_TriggerPadTracks_reco_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.centroid_",
                           TriggerPadTracks_reco_centroid_,
                           &b_TriggerPadTracks_reco_centroid_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.centroidX_",
                           TriggerPadTracks_reco_centroidX_,
                           &b_TriggerPadTracks_reco_centroidX_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.centroidY_",
                           TriggerPadTracks_reco_centroidY_,
                           &b_TriggerPadTracks_reco_centroidY_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.centroidZ_",
                           TriggerPadTracks_reco_centroidZ_,
                           &b_TriggerPadTracks_reco_centroidZ_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.residual_",
                           TriggerPadTracks_reco_residual_,
                           &b_TriggerPadTracks_reco_residual_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.residualX_",
                           TriggerPadTracks_reco_residualX_,
                           &b_TriggerPadTracks_reco_residualX_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.residualY_",
                           TriggerPadTracks_reco_residualY_,
                           &b_TriggerPadTracks_reco_residualY_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.residualZ_",
                           TriggerPadTracks_reco_residualZ_,
                           &b_TriggerPadTracks_reco_residualZ_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.nClusters_",
                           TriggerPadTracks_reco_nClusters_,
                           &b_TriggerPadTracks_reco_nClusters_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.beamEfrac_",
                           TriggerPadTracks_reco_beamEfrac_,
                           &b_TriggerPadTracks_reco_beamEfrac_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.px_",
                           TriggerPadTracks_reco_px_,
                           &b_TriggerPadTracks_reco_px_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.py_",
                           TriggerPadTracks_reco_py_,
                           &b_TriggerPadTracks_reco_py_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.pz_",
                           TriggerPadTracks_reco_pz_,
                           &b_TriggerPadTracks_reco_pz_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.x_", TriggerPadTracks_reco_x_,
                           &b_TriggerPadTracks_reco_x_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.y_", TriggerPadTracks_reco_y_,
                           &b_TriggerPadTracks_reco_y_);
  fChain->SetBranchAddress("TriggerPadTracks_reco.z_", TriggerPadTracks_reco_z_,
                           &b_TriggerPadTracks_reco_z_);
  fChain->SetBranchAddress("SiStripHits_reco", &SiStripHits_reco_,
                           &b_SiStripHits_reco_);
  fChain->SetBranchAddress("SiStripHits_reco.adcValues_",
                           SiStripHits_reco_adcValues_,
                           &b_SiStripHits_reco_adcValues_);
  fChain->SetBranchAddress("SiStripHits_reco.time_", SiStripHits_reco_time_,
                           &b_SiStripHits_reco_time_);
  fChain->SetBranchAddress("name_", &name_, &b_Trigger_reco_name_);
  fChain->SetBranchAddress("pass_", &pass_, &b_Trigger_reco_pass_);
  fChain->SetBranchAddress("variables_", &variables_,
                           &b_Trigger_reco_variables_);
  fChain->SetBranchAddress("FindableTracks_reco", &FindableTracks_reco_,
                           &b_FindableTracks_reco_);
  fChain->SetBranchAddress("FindableTracks_reco.particleTrackID_",
                           FindableTracks_reco_particleTrackID_,
                           &b_FindableTracks_reco_particleTrackID_);
  fChain->SetBranchAddress("FindableTracks_reco.is4sFindable_",
                           FindableTracks_reco_is4sFindable_,
                           &b_FindableTracks_reco_is4sFindable_);
  fChain->SetBranchAddress("FindableTracks_reco.is3s1aFindable_",
                           FindableTracks_reco_is3s1aFindable_,
                           &b_FindableTracks_reco_is3s1aFindable_);
  fChain->SetBranchAddress("FindableTracks_reco.is2s2aFindable_",
                           FindableTracks_reco_is2s2aFindable_,
                           &b_FindableTracks_reco_is2s2aFindable_);
  fChain->SetBranchAddress("FindableTracks_reco.is2aFindable_",
                           FindableTracks_reco_is2aFindable_,
                           &b_FindableTracks_reco_is2aFindable_);
  fChain->SetBranchAddress("FindableTracks_reco.is2sFindable_",
                           FindableTracks_reco_is2sFindable_,
                           &b_FindableTracks_reco_is2sFindable_);
  fChain->SetBranchAddress("FindableTracks_reco.is3sFindable_",
                           FindableTracks_reco_is3sFindable_,
                           &b_FindableTracks_reco_is3sFindable_);

  // loop over the input files

  for (auto inputFile : m_cfg.fileList) {
    // add file to the input chain
    fChain->Add(inputFile.c_str());
    ACTS_DEBUG("Adding File" << inputFile << "to tree " << m_cfg.treeName);
  }

  m_events = fChain->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

ActsExamples::RootLdmxSimHitReader::~RootLdmxSimHitReader() {
  /*
    delete m_layerID;
    delete m_moduleID;
    delete m_edep;
    delete m_time;
    delete m_px;
    delete m_py;
    delete m_pz;
    delete m_energy;
    delete m_x;
    delete m_y;
    delete m_z;
    delete m_pathLength;
    delete m_trackID;
    delete m_pdgID;
  */
}

std::string ActsExamples::RootLdmxSimHitReader::name() const {
  return "RootLdmxSimHitReader";
}

std::pair<size_t, size_t> ActsExamples::RootLdmxSimHitReader::availableEvents()
    const {
  return {0u, m_events};
}

ActsExamples::ProcessCode ActsExamples::RootLdmxSimHitReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read ldmx sim hits...");
  ACTS_DEBUG("The input TChain pointer=" << fChain);
  ACTS_DEBUG("The number of events=" << m_events);

  if (fChain && context.eventNumber < m_events) {
    // Lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);

    // Collection to be written

    // TODO::FIx
    std::vector<LdmxSpacePoint> mCollection;

    for (size_t ib = 0; ib < m_cfg.batchSize; ++ib) {
      // Read the correct entry: batch size * event_number + ib
      ACTS_DEBUG("Reading Entry... ");
      fChain->GetEntry(m_cfg.batchSize * context.eventNumber + ib);
      ACTS_VERBOSE("Reading entry: " << m_cfg.batchSize * context.eventNumber +
                                            ib);

      // std::cout<<"The size of sim hits is: " <<
      // TaggerSimHits_v12_<<std::endl; std::cout<<std::endl;

      int rotatePoints = 1;

      // Loop over the sim hits
      for (size_t idh = 0; idh < (size_t)TaggerSimHits_v12_; idh++) {
        // std::cout<<"Found sim hit "<< TaggerSimHits_v12_id_[idh] <<std::endl;
        // std::cout<<"Layer ID " <<TaggerSimHits_v12_layerID_[idh]<<std::endl;
        // std::cout<<"(x,y,z,t)=("<<TaggerSimHits_v12_x_[idh]<<","<<TaggerSimHits_v12_y_[idh]<<","<<TaggerSimHits_v12_z_[idh]<<","<<TaggerSimHits_v12_time_[idh]<<std::endl;

        // For the tracking example I'll rotate these points.
        // z -> x
        // x -> y
        // y -> z
        // t -> t

        if (rotatePoints == 1) {
          float pT = TaggerSimHits_v12_px_[idh] * TaggerSimHits_v12_px_[idh] +
                     TaggerSimHits_v12_py_[idh] * TaggerSimHits_v12_py_[idh];
          pT = std::sqrt(pT) * (1. / 1000);  // convert to GeV
          // ACTS_INFO("px is " << TaggerSimHits_v12_px_[idh] << ". and py is "
          //                    << TaggerSimHits_v12_py_[idh] << " and pz is "
          //                    << TaggerSimHits_v12_pz_[idh] << " and pT is "
          //                    << pT)

          LdmxSpacePoint ldmx_sp(
              TaggerSimHits_v12_z_[idh], TaggerSimHits_v12_x_[idh],
              TaggerSimHits_v12_y_[idh], TaggerSimHits_v12_time_[idh],
              TaggerSimHits_v12_layerID_[idh], TaggerSimHits_v12_id_[idh],
              TaggerSimHits_v12_trackID_[idh], pT,
              TaggerSimHits_v12_pz_[idh] / 1000.);
          mCollection.push_back(std::move(ldmx_sp));
        }

        // z->x
        // y->y
        // x->z

        // else if (rotatePoints == 2) {
        //   LdmxSpacePoint ldmx_sp(
        //       TaggerSimHits_v12_z_[idh], TaggerSimHits_v12_y_[idh],
        //       TaggerSimHits_v12_x_[idh], TaggerSimHits_v12_time_[idh],
        //       TaggerSimHits_v12_layerID_[idh], TaggerSimHits_v12_id_[idh],
        //       TaggerSimHits_v12_trackID_[idh]);
        //   mCollection.push_back(std::move(ldmx_sp));

        // }

        // else {
        //   LdmxSpacePoint ldmx_sp(
        //       TaggerSimHits_v12_x_[idh], TaggerSimHits_v12_y_[idh],
        //       TaggerSimHits_v12_z_[idh], TaggerSimHits_v12_time_[idh],
        //       TaggerSimHits_v12_layerID_[idh], TaggerSimHits_v12_id_[idh],
        //       TaggerSimHits_v12_trackID_[idh]);
        //   mCollection.push_back(std::move(ldmx_sp));
        // }

      }  // loop on the hits

    }  // end loop on the batch

    // Write the collection to the EventStore
    context.eventStore.add(m_cfg.outputCollection, std::move(mCollection));

    // std::cout<<"= = = Event Batch Done = = = "<<std::endl;
  }  // check if TChain exists and still events to process

  // Return a success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
