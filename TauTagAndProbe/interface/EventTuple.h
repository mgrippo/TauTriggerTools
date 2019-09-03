/*! Definition of a tuple with all event information that is required for the tau analysis.
*/

#pragma once

#include "TauTriggerTools/Common/interface/SmartTree.h"
#include "TauTriggerTools/Common/interface/TauIdResults.h"
#include <Math/VectorUtil.h>

#define TAU_ID(name, pattern, has_raw, wp_list) VAR(uint16_t, name) VAR(Float_t, name##raw)

#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define EVENT_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, npv) /* number of primary vertices */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, npu) /* number of in-time pu interactions added to the event */ \
    /* PF MET variables */ \
    VAR2(Float_t, met_pt, met_phi) /* pt and phi of the MET */ \
    /* Tag muon variables */ \
    VAR4(Float_t, muon_pt, muon_eta, muon_phi, muon_mass) /* 4-momentum of the muon */ \
    VAR(Int_t, muon_charge) /* muon charge */ \
    VAR(Float_t, muon_iso) /* muon pfRel isolation */ \
    VAR(Float_t, muon_mt) /* muon transverse mass */ \
    /* Basic tau variables */ \
    VAR(UInt_t, tau_sel) /* how tau was selected */ \
    VAR4(Float_t, tau_pt, tau_eta, tau_phi, tau_mass) /* 4-momentum of the tau */ \
    VAR(Int_t, tau_charge) /* tau charge */ \
    VAR(Int_t, lepton_gen_match) /* matching with leptons on the generator level (see Htautau Twiki for details):
                                    Electron = 1, Muon = 2, TauElectron = 3, TauMuon = 4, Tau = 5, NoMatch = 6 */\
    VAR(Int_t, lepton_gen_charge) /* charge of the matched gen lepton */ \
    VAR4(Float_t, lepton_gen_vis_pt, lepton_gen_vis_eta, \
                  lepton_gen_vis_phi, lepton_gen_vis_mass) /* visible 4-momentum of the matched gen lepton */ \
    /* Tau ID variables */ \
    VAR(Int_t, tau_decayMode) /* tau decay mode */ \
    VAR(Bool_t, tau_oldDecayModeFinding) /* tau passed the old decay mode finding requirements */ \
    TAU_IDS() \
    /* Tau transverse impact paramters.
       See cmssw/RecoTauTag/RecoTau/plugins/PFTauTransverseImpactParameters.cc for details */ \
    VAR(Float_t, tau_dxy) /* tau signed transverse impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    VAR(Float_t, tau_ip3d) /* tau signed 3D impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_ip3d_error) /* uncertainty of the 3D impact parameter measurement */ \
    VAR(Float_t, tau_dz) /* tau dz of the leadChargedHadrCand wrt to the primary vertex */ \
    VAR(Float_t, tau_dz_error) /* uncertainty of the tau dz measurement */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(tau_trigger, Event, EventTuple, EVENT_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(tau_trigger, EventTuple, EVENT_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef EVENT_DATA
#undef TAU_ID

namespace tau_trigger {

template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }
template<>
constexpr float DefaultFillValue<float>() { return -999.; }
template<>
constexpr int DefaultFillValue<int>() { return -999; }

} // namespace tau_tuple
