#include "TauTriggerTools/Common/interface/PatHelpers.h"


namespace tau_trigger {

double MuonIsolation(const pat::Muon& muon)
{
    const double pfIso = muon.pfIsolationR04().sumChargedHadronPt
                         + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt
                                    + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt);
    return pfIso / muon.polarP4().pt();
}

} // namespace tau_trigger
