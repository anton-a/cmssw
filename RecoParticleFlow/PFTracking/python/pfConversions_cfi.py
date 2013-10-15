import FWCore.ParameterSet.Config as cms

#
#  Author: M. Gouzevitch base on N. Marinelli work
#
pfConversions = cms.EDProducer("PFConversionProducer", 
    inputConversions = cms.InputTag("cleanedConversionsForPF", ""),
    inputPrimaryVertices = cms.InputTag("offlinePrimaryVertices")  
)


