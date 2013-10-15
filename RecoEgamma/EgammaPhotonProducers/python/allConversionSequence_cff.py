import FWCore.ParameterSet.Config as cms

#
#
# Tracker only conversion producer
from RecoEgamma.EgammaPhotonProducers.allConversions_cfi import *
# filtering for PFConversions
from RecoEgamma.EgammaPhotonProducers.conversionCleanerForPF_cfi import *
#allConversionSequence = cms.Sequence(allConversions)
# make a cleaned subset for PF
allConversionSequence = cms.Sequence(allConversions*cleanedConversionsForPF)


