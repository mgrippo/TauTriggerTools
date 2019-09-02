/*! Creates tuple for tau analysis.
*/

#include "Compression.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauTriggerTools/Common/interface/GenTruthTools.h"
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
        isMC(cfg.getParameter<bool>("isMC")),
        requireGenMatch(cfg.getParameter<bool>("requireGenMatch")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(consumes<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        signalMuon_token(consumes<pat::MuonRefVector>(cfg.getParameter<edm::InputTag>("signalMuon"))),
        taus_token(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
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

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));

        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());
        //summaryTuple().numberOfProcessedEvents++;

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

        const auto& PV = vertices->at(0);
        eventTuple().pv_x = static_cast<float>(PV.position().x());
        eventTuple().pv_y = static_cast<float>(PV.position().y());
        eventTuple().pv_z = static_cast<float>(PV.position().z());
        eventTuple().pv_chi2 = static_cast<float>(PV.chi2());
        eventTuple().pv_ndof = static_cast<float>(PV.ndof());


        edm::Handle<pat::MuonRefVector> signalMuonCollection;
        event.getByToken(signalMuon_token, signalMuonCollection);

        edm::Handle<pat::TauCollection> taus;
        event.getByToken(taus_token, taus);

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
        if(isMC)
            event.getByToken(genParticles_token, hGenParticles);

        // auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

        eventTuple.Fill();
    }

private:
    const bool isMC, requireGenMatch;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<pat::MuonRefVector> signalMuon_token;
    edm::EDGetTokenT<pat::TauCollection> taus_token;

    const TupleProducerData* data;
    EventTuple& eventTuple;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerTupleProducer = tau_trigger::TupleProducer;
DEFINE_FWK_MODULE(TauTriggerTupleProducer);
