import re
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from TauTriggerTools.Common.ProduceHelpers import *

options = VarParsing('analysis')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Text file with a list of the input root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Prefix to add to input file names.")
options.register('outputTupleFile', 'eventTuple.root', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Event tuple file.")
options.register('skipEvents', -1, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Number of events to skip")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of events to process.")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JSON file with lumi mask.")
options.register('period', 'Run2018', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Data taking period")
options.register('isMC', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Data or MC")
options.register('requireGenMatch', False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Require matching at the generator level")
options.register('useCustomHLT', False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Use custom HLT paths")
options.register('runDeepTau', True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Run DeepTau IDs")
options.parseArguments()

processName = "TagAndProbe"
process = cms.Process(processName)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.GlobalTag.globaltag = getGlobalTag(options.period, options.isMC)
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.TFileService = cms.Service('TFileService', fileName=cms.string(options.outputTupleFile))

if len(options.inputFileList) > 0:
    readFileList(process.source.fileNames, options.inputFileList, options.fileNamePrefix)
elif len(options.inputFiles) > 0:
    addFilesToList(process.source.fileNames, options.inputFiles, options.fileNamePrefix)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
if options.maxEvents > 0:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents > 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
if len(options.eventList) > 0:
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventList.split(','))
if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

# Update electron ID according recommendations from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
ele_era = {
    'Run2016': '2016-Legacy',
    'Run2017': '2017-Nov17ReReco',
}
ele_era.update(dict.fromkeys(['Run2018', 'Run2018ABC', 'Run2018D'], '2018-Prompt'))
setupEgammaPostRecoSeq(process, runVID=True, runEnergyCorrections=False, era=ele_era[options.period])

# Update tau IDs according recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
updatedTauName = "slimmedTausNewID"
tauIdsToKeep = [ "2017v2" ]
if options.runDeepTau:
    tauIdsToKeep.append("deepTau2017v2p1")
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug=False, updatedTauName=updatedTauName,
                                          toKeep=tauIdsToKeep)
tauIdEmbedder.runTauID()
tauSrc_InputTag = cms.InputTag(updatedTauName)

# Update MET filters according recommendations from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# Using post-Moriond2019 (a more complete) list of noisy crystals
process.metFilterSequence = cms.Sequence()
customMetFilters = cms.PSet()
if options.period in [ 'Run2017', 'Run2018', 'Run2018ABC', 'Run2018D' ]:
    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
    baddetEcallist = cms.vuint32([
        872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336,
        872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339,
        872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185,
        872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694,
        872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609,
        872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413,
        872421414,872421031,872423083,872421439])
    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt = cms.double(50.),
        baddetEcal = baddetEcallist,
        taggingMode = cms.bool(True),
        debug = cms.bool(False)
    )
    process.metFilterSequence += process.ecalBadCalibReducedMINIAODFilter
    customMetFilters.ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")

# Re-apply MET corrections
if options.period in [ 'Run2016', 'Run2017' ]:
    met_run_params = { }
    if options.period == 'Run2017':
        met_run_params = {
            'fixEE2017': True,
            'fixEE2017Params': {
                'userawPt': True,
                'ptThreshold':50.0,
                'minEtaThreshold':2.65,
                'maxEtaThreshold': 3.139
            }
        }
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process, isData = not options.isMC, **met_run_params)
    metInputTag = cms.InputTag('slimmedMETs', '', processName)
else:
    metInputTag = cms.InputTag('slimmedMETs')

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = ['HLT_IsoMu27_v*'],
    andOr = cms.bool(True), # how to deal with multiple triggers:
                            # True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid
)

process.selectionFilter = cms.EDFilter("TauTriggerSelectionFilter",
    electrons        = cms.InputTag('slimmedElectrons'),
    muons            = cms.InputTag('slimmedMuons'),
    jets             = cms.InputTag('slimmedJets'),
    met              = metInputTag,
    triggerResults   = cms.InputTag('TriggerResults', '', 'PAT'),
    customMetFilters = customMetFilters,
    btagThreshold    = cms.double(-1),
    metFilters       = cms.vstring(getMetFilters(options.period, options.isMC)),
    mtCut            = cms.double(-1)
)

process.tupleProducer = cms.EDProducer("TauTriggerTupleProducer",
    isMC            = cms.bool(options.isMC),
    requireGenMatch = cms.bool(options.requireGenMatch),
    genEvent        = cms.InputTag('generator'),
    puInfo          = cms.InputTag('slimmedAddPileupInfo'),
    genParticles    = cms.InputTag('prunedGenParticles'),
    vertices        = cms.InputTag('offlineSlimmedPrimaryVertices'),
    signalMuon      = cms.InputTag('selectionFilter'),
    taus            = tauSrc_InputTag,
    jets            = cms.InputTag('slimmedJets'),
    met             = metInputTag,
    btagThreshold   = cms.double(getBtagThreshold(options.period, 'Loose')),
)

process.p = cms.Path(
    process.egammaPostRecoSeq +
    process.metFilterSequence +
    process.hltFilter +
    process.selectionFilter +
    process.rerunMvaIsolationSequence +
    getattr(process, updatedTauName) +
    process.tupleProducer
)

# Verbosity customization
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = getReportInterval(process.maxEvents.input.value())
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
