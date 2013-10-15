import FWCore.ParameterSet.Config as cms

# filter allConversions for use in PF

cleanedConversionsForPF = cms.EDProducer("ConversionCleanerForPF", 
    inputConversionsCollection = cms.InputTag("allConversions", "")
#replace some of the hardcoded values with parameters
)
