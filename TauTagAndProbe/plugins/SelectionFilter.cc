/*! Apply tau trigger selection vetoes.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TauTriggerTools/Common/interface/AnalysisTypes.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"

namespace tau_trigger {

class SelectionFilter : public edm::EDFilter {
public:
    SelectionFilter(const edm::ParameterSet& cfg) :
        btagThreshold(cfg.getParameter<double>("btagThreshold")),
        mtCut(cfg.getParameter<double>("mtCut")),
        metFilters(cfg.getParameter<std::vector<std::string>>("metFilters")),
        electrons_token(consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("electrons"))),
        muons_token(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        jets_token(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
        met_token(consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"))),
        triggerResults_token(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults")))
    {
        const edm::ParameterSet& customMetFilters = cfg.getParameterSet("customMetFilters");
        for(const auto& filterName : customMetFilters.getParameterNames()) {
            customMetFilters_token[filterName] =
                mayConsume<bool>(customMetFilters.getParameter<edm::InputTag>(filterName));
        }
        produces<pat::MuonRefVector>();
    }

private:
    virtual bool filter(edm::Event& event, const edm::EventSetup&) override
    {
        edm::Handle<pat::MuonCollection> muons;
        event.getByToken(muons_token, muons);

        // Find signal muon
        std::vector<pat::MuonRef> signalMuonCandidates;
        for(size_t n = 0; n < muons->size(); ++n) {
            const pat::Muon& muon = muons->at(n);
            if(muon.polarP4().pt() > 2.4 && std::abs(muon.polarP4().eta()) < 2.1 && muon.isMediumMuon())
                signalMuonCandidates.emplace_back(muons, n);
        }
        if(signalMuonCandidates.empty()) return false;
        std::sort(signalMuonCandidates.begin(), signalMuonCandidates.end(),
                  [](const pat::MuonRef& a, const pat::MuonRef& b) { return MuonIsolation(*a) < MuonIsolation(*b); });
        const pat::Muon& signalMuon = *signalMuonCandidates.at(0);

        // Apply third lepton veto
        for(const pat::Muon& muon : *muons) {
            if(&muon != &signalMuon && muon.isLooseMuon() && muon.polarP4().pt() > 10
                    && std::abs(muon.polarP4().eta()) < 2.4 && MuonIsolation(muon) < 0.3)
                return false;
        }
        edm::Handle<pat::ElectronCollection> electrons;
        event.getByToken(electrons_token, electrons);
        for(const pat::Electron& ele : *electrons) {
            if(ele.polarP4().pt() > 10 && std::abs(ele.polarP4().eta()) < 2.5
                    && ele.electronID("mvaEleID-Fall17-iso-V2-wpLoose") > 0.5)
                return false;
        }

        // Apply MT cut (if enabled)
        if(mtCut > 0) {
            edm::Handle<pat::METCollection> metCollection;
            event.getByToken(met_token, metCollection);
            const pat::MET& met = metCollection->at(0);
            const analysis::LorentzVectorM met_p4(met.pt(), 0, met.phi(), 0);
            if(Calculate_MT(signalMuon.polarP4(), met_p4) > mtCut) return false;
        }

        // Apply b tag veto (if enabled)
        if(btagThreshold > 0) {
            edm::Handle<pat::JetCollection> jets;
            event.getByToken(jets_token, jets);
            for(const pat::Jet& jet : *jets) {
                const auto btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")
                                  + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")
                                  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
                if(jet.polarP4().pt() > 20 && std::abs(jet.polarP4().eta()) < 2.4 && btag > btagThreshold)
                    return false;
            }
        }

        // Apply MET filters
        edm::Handle<edm::TriggerResults> triggerResults;
        event.getByToken(triggerResults_token, triggerResults);
        const edm::TriggerNames& triggerNames = event.triggerNames(*triggerResults);
        const auto passFilter = [&](const std::string& filter) {
            auto iter = customMetFilters_token.find(filter);
            if(iter != customMetFilters_token.end()) {
                edm::Handle<bool> result;
                event.getByToken(iter->second, result);
                return *result;
            }
            const size_t index = triggerNames.triggerIndex(filter);
            if(index == triggerNames.size())
                throw cms::Exception("TauTriggerSelectionFilter") << "MET filter '" << filter << "' not found.";
            return triggerResults->accept(index);
        };
        for(const std::string& metFilter : metFilters) {
            if(!passFilter(metFilter))
                return false;
        }

        // Put the signal muon into the event
        auto signalMuonOutput = std::make_unique<pat::MuonRefVector>();
        signalMuonOutput->push_back(signalMuonCandidates.at(0));
        event.put(std::move(signalMuonOutput));
        return true;
    }

private:
    const double btagThreshold, mtCut;
    const std::vector<std::string> metFilters;

    edm::EDGetTokenT<pat::ElectronCollection> electrons_token;
    edm::EDGetTokenT<pat::MuonCollection> muons_token;
    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::METCollection> met_token;
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_token;
    std::map<std::string, edm::EDGetTokenT<bool>> customMetFilters_token;

};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerSelectionFilter = tau_trigger::SelectionFilter;
DEFINE_FWK_MODULE(TauTriggerSelectionFilter);
