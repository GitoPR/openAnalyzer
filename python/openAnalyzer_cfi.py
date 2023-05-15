import FWCore.ParameterSet.Config as cms

openAnalyzer = cms.EDAnalyzer("openAnalyzer",
                               ecalhits = cms.InputTag("particleFlowRecHitECAL", "Cleaned", "RECO")
)
