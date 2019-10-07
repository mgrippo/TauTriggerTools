/*! Apply tau trigger selection vetoes.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TauTriggerTools/Common/interface/AnalysisTypes.h"
#include "TauTriggerTools/Common/interface/CutTools.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "TauTriggerTools/Common/interface/GenTruthTools.h"

class GenTauFilter : public edm::EDFilter {
public:

    GenTauFilter(const edm::ParameterSet& cfg) :
        genParticles_token(consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        n_tau_had(cfg.getParameter<int>("n_tau_had")),
        n_tau_e(cfg.getParameter<int>("n_tau_e")),
        n_tau_mu(cfg.getParameter<int>("n_tau_mu"))
    { }

private:
    virtual bool filter(edm::Event& event, const edm::EventSetup&) override
    {
        bool result = true;
        edm::Handle<std::vector<reco::GenParticle>> genParticles;
        event.getByToken(genParticles_token, genParticles);
        std::vector<analysis::gen_truth::LeptonMatchResult> lepton_results = analysis::gen_truth::CollectGenLeptons(*genParticles);
        int count_tau_h = 0;
        int count_tau_e = 0;
        int count_tau_mu = 0;
        for(unsigned n = 0; n < lepton_results.size(); ++n){
            analysis::gen_truth::LeptonMatchResult lepton_result = lepton_results.at(n);
            if(lepton_result.match == analysis::GenLeptonMatch::Tau) ++count_tau_h;
            if(lepton_result.match == analysis::GenLeptonMatch::TauElectron) ++count_tau_e;
            if(lepton_result.match == analysis::GenLeptonMatch::TauMuon) ++count_tau_mu;
        }
        //std::cout << "tau_h: " << count_tau_h << ", tau_e: " << count_tau_e << ", tau_mu: " << count_tau_mu << std::endl;
        if(count_tau_h != n_tau_had || count_tau_e != n_tau_e || count_tau_mu != n_tau_mu)
            result = false;
        return result;
    }

    void endJob()
    {}

private:
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    int n_tau_had,n_tau_e,n_tau_mu;

};


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenTauFilter);
