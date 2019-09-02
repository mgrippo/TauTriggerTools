#pragma once

#include <Math/VectorUtil.h>
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace tau_trigger {

double MuonIsolation(const pat::Muon& muon);

template<typename LVector1, typename LVector2>
double Calculate_MT(const LVector1& lepton_p4, const LVector2& met_p4)
{
    const double delta_phi = ROOT::Math::VectorUtil::DeltaPhi(lepton_p4, met_p4);
    return std::sqrt( 2.0 * lepton_p4.Pt() * met_p4.Pt() * ( 1.0 - std::cos(delta_phi) ) );
}

}
