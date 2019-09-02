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
    VAR(Float_t, rho) /* fixed grid energy density */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, npu) /* number of in-time pu interactions added to the event */ \
    VAR3(Float_t, pv_x, pv_y, pv_z) /* position of the primary vertex (PV) */ \
    VAR(Float_t, pv_chi2) /* chi^2 of the primary vertex (PV) */ \
    VAR(Float_t, pv_ndof) /* number of degrees of freedom of the primary vertex (PV) */ \
    /* Basic tau variables */ \
    VAR(Int_t, tau_sel) /* how tau was selected */ \
    VAR4(Float_t, tau_pt, tau_eta, tau_phi, tau_mass) /* 4-momentum of the tau */ \
    VAR(Int_t, tau_charge) /* tau charge */ \
    VAR(Int_t, lepton_gen_match) /* matching with leptons on the generator level (see Htautau Twiki for details):
                                    Electron = 1, Muon = 2, TauElectron = 3, TauMuon = 4, Tau = 5, NoMatch = 6 */\
    VAR(Int_t, lepton_gen_charge) /* charge of the matched gen lepton */ \
    VAR4(Float_t, lepton_gen_pt, lepton_gen_eta, \
                  lepton_gen_phi, lepton_gen_mass) /* 4-momentum of the matched gen lepton */ \
    VAR(std::vector<Int_t>, lepton_gen_vis_pdg) /* PDG of the matched lepton */ \
    VAR4(std::vector<Float_t>, lepton_gen_vis_pt, lepton_gen_vis_eta, \
                               lepton_gen_vis_phi, lepton_gen_vis_mass) /* 4-momenta of the visible products
                                                                           of the matched gen lepton */ \
    VAR(Int_t, qcd_gen_match) /* matching with QCD particles on the generator level:
                                 NoMatch = 0, Down = 1, Up = 2, Strange = 3, Charm = 4, Bottom = 5, Top = 6,
                                 Gluon = 21 */ \
    VAR(Int_t, qcd_gen_charge) /* charge of the matched gen QCD particle */ \
    VAR4(Float_t, qcd_gen_pt, qcd_gen_eta, qcd_gen_phi, qcd_gen_mass) /* 4-momentum of the matched gen QCD particle */ \
    /* Tau ID variables */ \
    VAR(Int_t, tau_decayMode) /* tau decay mode */ \
    VAR(Int_t, tau_decayModeFinding) /* tau passed the old decay mode finding requirements */ \
    VAR(Int_t, tau_decayModeFindingNewDMs) /* tau passed the new decay mode finding requirements */ \
    VAR(Float_t, chargedIsoPtSum) /* sum of the transverse momentums of charged pf candidates inside
                                     the tau isolation cone with dR < 0.5 */ \
    VAR(Float_t, chargedIsoPtSumdR03) /* sum of the transverse momentums of charged pf candidates inside
                                         the tau isolation cone with dR < 0.3 */ \
    VAR(Float_t, footprintCorrection) /* tau footprint correction inside the tau isolation cone with dR < 0.5 */ \
    VAR(Float_t, footprintCorrectiondR03) /* tau footprint correction inside the tau isolation cone with dR < 0.3 */ \
    VAR(Float_t, neutralIsoPtSum) /* sum of the transverse momentums of neutral pf candidates inside
                                     the tau isolation cone with dR < 0.5 */ \
    VAR(Float_t, neutralIsoPtSumWeight) /* weighted sum of the transverse momentums of neutral pf candidates inside
                                           the tau isolation cone with dR < 0.5 */ \
    VAR(Float_t, neutralIsoPtSumWeightdR03) /* weighted sum of the transverse momentums of neutral pf candidates inside
                                               the tau isolation cone with dR < 0.3 */ \
    VAR(Float_t, neutralIsoPtSumdR03) /* sum of the transverse momentums of neutral pf candidates inside
                                         the tau isolation cone with dR < 0.3 */ \
    VAR(Float_t, photonPtSumOutsideSignalCone) /* sum of the transverse momentums of photons
                                                  inside the tau isolation cone with dR < 0.5 */ \
    VAR(Float_t, photonPtSumOutsideSignalConedR03) /* sum of the transverse momentums of photons inside
                                                      the tau isolation cone with dR < 0.3 */ \
    VAR(Float_t, puCorrPtSum) /* pile-up correction for the sum of the transverse momentums */ \
    TAU_IDS() \
    /* Tau transverse impact paramters.
       See cmssw/RecoTauTag/RecoTau/plugins/PFTauTransverseImpactParameters.cc for details */ \
    VAR3(Float_t, tau_dxy_pca_x, tau_dxy_pca_y, tau_dxy_pca_z) /* The point of closest approach (PCA) of
                                                                  the leadPFChargedHadrCand to the primary vertex */ \
    VAR(Float_t, tau_dxy) /* tau signed transverse impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    VAR(Float_t, tau_ip3d) /* tau signed 3D impact parameter wrt to the primary vertex */ \
    VAR(Float_t, tau_ip3d_error) /* uncertainty of the 3D impact parameter measurement */ \
    VAR(Float_t, tau_dz) /* tau dz of the leadChargedHadrCand wrt to the primary vertex */ \
    VAR(Float_t, tau_dz_error) /* uncertainty of the tau dz measurement */ \
    VAR(Int_t, tau_hasSecondaryVertex) /* tau has the secondary vertex */ \
    VAR3(Float_t, tau_sv_x, tau_sv_y, tau_sv_z) /* position of the secondary vertex */ \
    VAR3(Float_t, tau_flightLength_x, tau_flightLength_y, tau_flightLength_z) /* flight length of the tau */ \
    VAR(Float_t, tau_flightLength_sig) /* significance of the flight length measurement */ \
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
