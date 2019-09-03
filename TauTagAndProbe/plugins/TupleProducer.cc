/*! Creates tuple for tau analysis.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "Compression.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "TauTriggerTools/TauTagAndProbe/interface/EventTuple.h"

namespace tau_trigger {

struct TupleProducerData {
    using Mutex = EventTuple::Mutex;
    using LockGuard = std::lock_guard<Mutex>;

    std::unique_ptr<EventTuple> eventTuple;

public:
    TupleProducerData(TFile& file)
    {
        eventTuple = std::make_unique<EventTuple>("events", &file, false);
    }
};

class TupleProducer : public edm::stream::EDProducer<edm::GlobalCache<TupleProducerData>> {
public:
    TupleProducer(const edm::ParameterSet& cfg, const TupleProducerData* producerData) :
        btagThreshold(cfg.getParameter<double>("btagThreshold")),
        isMC(cfg.getParameter<bool>("isMC")),
        requireGenMatch(cfg.getParameter<bool>("requireGenMatch")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(consumes<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        signalMuon_token(consumes<pat::MuonRefVector>(cfg.getParameter<edm::InputTag>("signalMuon"))),
        taus_token(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        jets_token(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
        met_token(consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"))),
        data(producerData),
        eventTuple(*data->eventTuple)
    {
        produces<bool>();
    }

    static std::unique_ptr<TupleProducerData> initializeGlobalCache(const edm::ParameterSet&)
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kLZ4);
        file.SetCompressionLevel(4);
        return std::make_unique<TupleProducerData>(file);
    }

    static void globalEndJob(TupleProducerData* data)
    {
        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());
        data->eventTuple->Write();
    }

private:
    static constexpr float default_value = ::tau_trigger::DefaultFillValue<float>();
    static constexpr int default_int_value = ::tau_trigger::DefaultFillValue<int>();
    static constexpr double deltaR2Thr = 0.5*0.5;

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));

        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());

        eventTuple().run  = event.id().run();
        eventTuple().lumi = event.id().luminosityBlock();
        eventTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::Vertex>> vertices;
        event.getByToken(vertices_token, vertices);
        eventTuple().npv = static_cast<int>(vertices->size());

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            eventTuple().genEventWeight = static_cast<float>(genEvent->weight());

            edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
            event.getByToken(puInfo_token, puInfo);
            eventTuple().npu = ::analysis::gen_truth::GetNumberOfPileUpInteractions(puInfo);
        }

        edm::Handle<pat::MuonRefVector> signalMuonCollection;
        event.getByToken(signalMuon_token, signalMuonCollection);
        const pat::Muon* muon = signalMuonCollection.isValid() && !signalMuonCollection->empty()
                              ? &(*signalMuonCollection->at(0)) : nullptr;

        edm::Handle<pat::METCollection> metCollection;
        event.getByToken(met_token, metCollection);
        const pat::MET& met = metCollection->at(0);
        const analysis::LorentzVectorM met_p4(met.pt(), 0, met.phi(), 0);
        eventTuple().met_pt = static_cast<float>(met.pt());
        eventTuple().met_phi = static_cast<float>(met.phi());
        eventTuple().muon_pt = muon ? static_cast<float>(muon->polarP4().pt()) : default_value;
        eventTuple().muon_eta = muon ? static_cast<float>(muon->polarP4().eta()) : default_value;
        eventTuple().muon_phi = muon ? static_cast<float>(muon->polarP4().phi()) : default_value;
        eventTuple().muon_mass = muon ? static_cast<float>(muon->polarP4().mass()) : default_value;
        eventTuple().muon_charge = muon ? muon->charge() : default_int_value;
        eventTuple().muon_iso = muon ? MuonIsolation(*muon) : default_value;
        eventTuple().muon_mt = muon ? Calculate_MT(muon->polarP4(), analysis::LorentzVectorM(met.pt(), 0, met.phi(), 0))
                                    : default_value;

        edm::Handle<pat::TauCollection> taus;
        event.getByToken(taus_token, taus);

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
        if(isMC)
            event.getByToken(genParticles_token, hGenParticles);
        auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

        edm::Handle<pat::JetCollection> jets;
        event.getByToken(jets_token, jets);

        const auto& selected_taus = CollectTaus(muon, *taus);
        for(const auto& tau_entry : selected_taus) {
            const pat::Tau& tau = *tau_entry.first;
            const unsigned selection = tau_entry.second;
            if(!PassBtagVeto(muon, tau, *jets)) continue;

            eventTuple().tau_sel = selection;
            eventTuple().tau_pt = static_cast<float>(tau.polarP4().pt());
            eventTuple().tau_eta = static_cast<float>(tau.polarP4().eta());
            eventTuple().tau_phi = static_cast<float>(tau.polarP4().phi());
            eventTuple().tau_mass = static_cast<float>(tau.polarP4().mass());
            eventTuple().tau_charge = tau.charge();

            analysis::gen_truth::LeptonMatchResult leptonMatch;
            if(genParticles)
                leptonMatch = analysis::gen_truth::LeptonGenMatch(tau.polarP4(), *genParticles);
            if(requireGenMatch && leptonMatch.match != analysis::GenLeptonMatch::Tau) continue;
            const bool has_lepton = leptonMatch.match != analysis::GenLeptonMatch::NoMatch;
            eventTuple().lepton_gen_match = static_cast<int>(leptonMatch.match);
            eventTuple().lepton_gen_charge = has_lepton ? leptonMatch.gen_particle->charge() : default_int_value;
            eventTuple().lepton_gen_vis_pt = has_lepton ? static_cast<float>(leptonMatch.visible_daughters_p4.pt())
                                                        : default_value;
            eventTuple().lepton_gen_vis_eta = has_lepton ? static_cast<float>(leptonMatch.visible_daughters_p4.eta())
                                                         : default_value;
            eventTuple().lepton_gen_vis_phi = has_lepton ? static_cast<float>(leptonMatch.visible_daughters_p4.phi())
                                                         : default_value;
            eventTuple().lepton_gen_vis_mass = has_lepton ? static_cast<float>(leptonMatch.visible_daughters_p4.mass())
                                                          : default_value;

            eventTuple().tau_decayMode = tau.decayMode();
            eventTuple().tau_oldDecayModeFinding = tau.tauID("decayModeFinding") > 0.5f;

            for(const auto& tau_id_entry : analysis::tau_id::GetTauIdDescriptors()) {
                const auto& desc = tau_id_entry.second;
                desc.FillTuple(eventTuple, &tau, default_value);
            }

            eventTuple().tau_dxy = tau.dxy();
            eventTuple().tau_dxy_error = tau.dxy_error();
            eventTuple().tau_ip3d = tau.ip3d();
            eventTuple().tau_ip3d_error = tau.ip3d_error();

            auto leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
            eventTuple().tau_dz = leadChargedHadrCand ? leadChargedHadrCand->dz() : default_value;
            eventTuple().tau_dz_error = leadChargedHadrCand && leadChargedHadrCand->hasTrackDetails()
                    ? leadChargedHadrCand->dzError() : default_value;

            eventTuple.Fill();
        }
    }

    std::map<const pat::Tau*, unsigned> CollectTaus(const pat::Muon* muon, const pat::TauCollection& taus) const
    {
        std::map<const pat::Tau*, unsigned> selected_taus;
        std::vector<const pat::Tau*> tau_vec;
        for(const auto& tau : taus) {
            auto leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
            if(tau.polarP4().pt() > 20 && std::abs(tau.polarP4().eta()) < 2.3
                    && leadChargedHadrCand && std::abs(leadChargedHadrCand->dz()) < 0.2
                    && (!muon || reco::deltaR2(muon->polarP4(), tau.polarP4()) > deltaR2Thr)) {
                tau_vec.push_back(&tau);
            }
        }
        if(!tau_vec.empty()) {
            const pat::Tau* tau = *std::max_element(tau_vec.begin(), tau_vec.end(),
                [](const pat::Tau* tau1, const pat::Tau* tau2) { return tau1->polarP4().pt() < tau2->polarP4().pt(); });
            selected_taus[tau] |= static_cast<unsigned>(analysis::TauSelection::pt);

            const auto selectBest = [&](const std::string& discr_name, analysis::TauSelection sel) {
                if(tau_vec.at(0)->isTauIDAvailable(discr_name)) {
                    tau = *std::max_element(tau_vec.begin(), tau_vec.end(),
                        [&](const pat::Tau* tau1, const pat::Tau* tau2) {
                            return tau1->tauID(discr_name) < tau2->tauID(discr_name);
                        });
                    selected_taus[tau] |= static_cast<unsigned>(sel);
                }
            };

            selectBest("byIsolationMVArun2017v2DBoldDMwLTraw2017", analysis::TauSelection::MVA);
            selectBest("byDeepTau2017v2p1VSjetraw", analysis::TauSelection::DeepTau);
        }
        return selected_taus;
    }

    bool PassBtagVeto(const pat::Muon* muon, const pat::Tau& tau, const pat::JetCollection& jets) const
    {
        if(btagThreshold > 0) {
            for(const pat::Jet& jet : jets) {
                const auto btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")
                                  + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")
                                  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
                if((!muon || reco::deltaR2(muon->polarP4(), jet.polarP4()) > deltaR2Thr)
                        && reco::deltaR2(tau.polarP4(), jet.polarP4()) > deltaR2Thr
                        && jet.polarP4().pt() > 20 && std::abs(jet.polarP4().eta()) < 2.4
                        && btag > btagThreshold)
                    return false;
            }
        }
        return true;
    }

private:
    const double btagThreshold;
    const bool isMC, requireGenMatch;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<pat::MuonRefVector> signalMuon_token;
    edm::EDGetTokenT<pat::TauCollection> taus_token;
    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::METCollection> met_token;

    const TupleProducerData* data;
    EventTuple& eventTuple;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerTupleProducer = tau_trigger::TupleProducer;
DEFINE_FWK_MODULE(TauTriggerTupleProducer);
