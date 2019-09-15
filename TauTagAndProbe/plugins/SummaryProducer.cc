/*! Creates tuple for tau analysis.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "Compression.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Common/interface/TriggerDescriptor.h"
#include "TauTriggerTools/TauTagAndProbe/interface/SummaryTuple.h"

namespace tau_trigger {

using namespace analysis;

struct SummaryProducerData {
    using Mutex = SummaryTuple::Mutex;
    using LockGuard = std::lock_guard<Mutex>;
    using clock = std::chrono::system_clock;

    const clock::time_point start;
    std::unique_ptr<SummaryTuple> summaryTuple;
    std::unique_ptr<ExpressTuple> expressTuple;

    SummaryProducerData(TFile& file, bool createExpressTuple) :
        start(clock::now()),
        summaryTuple(std::make_unique<SummaryTuple>("summary", &file, false))
    {
        if(createExpressTuple)
            expressTuple = std::make_unique<ExpressTuple>("all_events", &file, false);
    }
};

class SummaryProducer : public edm::stream::EDProducer<edm::GlobalCache<SummaryProducerData>> {
public:
    SummaryProducer(const edm::ParameterSet& cfg, const SummaryProducerData* data) :
        isMC(cfg.getParameter<bool>("isMC")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(mayConsume<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        summaryTuple(*data->summaryTuple),
        expressTuple(data->expressTuple.get())
    {
        produces<bool>();
    }

    static std::unique_ptr<SummaryProducerData> initializeGlobalCache(const edm::ParameterSet& cfg)
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kLZ4);
        file.SetCompressionLevel(4);
        const bool isMC = cfg.getParameter<bool>("isMC");
        auto data = std::make_unique<SummaryProducerData>(file, isMC);
        (*data->summaryTuple)().numberOfProcessedEvents = 0;
        (*data->summaryTuple)().totalGenEventWeight = 0;

        TriggerDescriptorCollection hltPaths(cfg.getParameter<edm::VParameterSet>("hltPaths"));
        for(unsigned n = 0; n < hltPaths.size(); ++n) {
            (*data->summaryTuple)().trigger_index.push_back(n);
            (*data->summaryTuple)().trigger_pattern.push_back(hltPaths.at(n).path);
        }
        return data;
    }

    static void globalEndJob(SummaryProducerData* data)
    {
        SummaryProducerData::LockGuard lock(data->summaryTuple->GetMutex());
        if(data->expressTuple)
            data->expressTuple->Write();
        const auto stop = SummaryProducerData::clock::now();
        (*data->summaryTuple)().exeTime = std::chrono::duration_cast<std::chrono::seconds>(stop - data->start).count();
        data->summaryTuple->Fill();
        data->summaryTuple->Write();
    }

private:
    static constexpr float default_value = ::tau_trigger::DefaultFillValue<float>();
    static constexpr int default_int_value = ::tau_trigger::DefaultFillValue<int>();

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));

        SummaryProducerData::LockGuard lock(summaryTuple.GetMutex());

        summaryTuple().numberOfProcessedEvents++;

        if(expressTuple) {
            (*expressTuple)().run  = event.id().run();
            (*expressTuple)().lumi = event.id().luminosityBlock();
            (*expressTuple)().evt  = event.id().event();

            edm::Handle<std::vector<reco::Vertex>> vertices;
            event.getByToken(vertices_token, vertices);
            (*expressTuple)().npv = static_cast<int>(vertices->size());
        }

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            const float genWeight = static_cast<float>(genEvent->weight());
            summaryTuple().totalGenEventWeight += genWeight;
            if(expressTuple) {
                (*expressTuple)().genEventWeight = genWeight;

                edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
                event.getByToken(puInfo_token, puInfo);
                (*expressTuple)().npu = gen_truth::GetNumberOfPileUpInteractions(puInfo);
            }
        } else if(expressTuple) {
            (*expressTuple)().genEventWeight = default_int_value;
            (*expressTuple)().npu = default_int_value;
        }

        if(expressTuple)
            expressTuple->Fill();
    }

private:
    const bool isMC;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;

    SummaryTuple& summaryTuple;
    ExpressTuple* expressTuple;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerSummaryTupleProducer = tau_trigger::SummaryProducer;
DEFINE_FWK_MODULE(TauTriggerSummaryTupleProducer);
