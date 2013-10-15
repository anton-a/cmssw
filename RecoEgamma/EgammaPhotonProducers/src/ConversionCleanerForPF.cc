// -*- C++ -*-
//
// Package:    EgammaPhotonProducers
// Class:      ConversionCleanerForPF
// 
// Created:  Mon Aug 26 15:22:31 CET 2013
//           Anton Anastassov
//
// NOTE: The first version contains exactly the same filtering code as implemented in the 
//       pre-CMSSW_7XX PFConversionProducer

#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionCleanerForPF.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionLikelihoodCalculator.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionHitChecker.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

using namespace std;
using namespace edm;

ConversionCleanerForPF::ConversionCleanerForPF(const ParameterSet& iConfig) {

  produces<reco::ConversionRefVector>();
  inputConversions_ = 
    consumes<reco::ConversionCollection>(iConfig.getParameter<InputTag>("inputConversionsCollection"));
}

ConversionCleanerForPF::~ConversionCleanerForPF()
{
}

void ConversionCleanerForPF::produce(Event& iEvent, const EventSetup& iSetup)
{
  
  //create the empty collection 
  auto_ptr<reco::ConversionRefVector> cleanedConversions (new reco::ConversionRefVector);

  Handle<reco::ConversionCollection> convCollHandle;
  iEvent.getByToken(inputConversions_, convCollHandle);
  const reco::ConversionCollection& convColl = *(convCollHandle.product());
  
   // (AA) looks like the code below can be optimized and cleaned-up 
   // consult with the original authors
   
  // CLEAN CONVERSION COLLECTION FOR DUPLICATES     
  for( unsigned int icoll1=0; icoll1 < convColl.size(); icoll1++) 
    { 
      if (( !convColl[icoll1].quality(reco::Conversion::arbitratedMergedEcalGeneral)) || (!convColl[icoll1].quality(reco::Conversion::highPurity))) continue;
      
      bool greater_prob=false;
      std::vector<edm::RefToBase<reco::Track> > tracksRefColl1 = convColl[icoll1].tracks();      
      for(unsigned it1 = 0; it1 < tracksRefColl1.size(); it1++)
  {
    reco::TrackRef trackRef1 = (tracksRefColl1[it1]).castTo<reco::TrackRef>();
   
    for( unsigned int icoll2=0; icoll2 < convColl.size(); icoll2++) 
      {
        if(icoll1==icoll2)continue;
        if (( !convColl[icoll2].quality(reco::Conversion::arbitratedMergedEcalGeneral)) || (!convColl[icoll2].quality(reco::Conversion::highPurity))) continue;
        std::vector<edm::RefToBase<reco::Track> > tracksRefColl2 = convColl[icoll2].tracks();     
        for(unsigned it2 = 0; it2 < tracksRefColl2.size(); it2++)
    {
      reco::TrackRef trackRef2 = (tracksRefColl2[it2]).castTo<reco::TrackRef>();
      double like1=-999;
      double like2=-999;
      //number of shared hits
      int shared=0;
      for(trackingRecHit_iterator iHit1=trackRef1->recHitsBegin(); iHit1!=trackRef1->recHitsEnd(); iHit1++) 
        {
          const TrackingRecHit *h_1=iHit1->get();
          if(h_1->isValid()){     
      for(trackingRecHit_iterator iHit2=trackRef2->recHitsBegin(); iHit2!=trackRef2->recHitsEnd(); iHit2++)
        {
          const TrackingRecHit *h_2=iHit2->get();
          if(h_2->isValid() && h_1->sharesInput(h_2, TrackingRecHit::some))shared++;//count number of shared hits
        }
          }
        }     
      float frac=0;
      //number of valid hits in tracks that are duplicates
      float size1=trackRef1->found();
      float size2=trackRef2->found();
      //divide number of shared hits by the total number of hits for the track with less hits
      if(size1>size2)frac=(double)shared/size2;
      else frac=(double)shared/size1;
      if(frac>0.9)
        {
          like1=ChiSquaredProbability(convColl[icoll1].conversionVertex().chi2(), convColl[icoll1].conversionVertex().ndof());
          like2=ChiSquaredProbability(convColl[icoll2].conversionVertex().chi2(), convColl[icoll2].conversionVertex().ndof());
        }
      if(like2>like1)
        {greater_prob=true;  break;}
    }//end loop over tracks in collection 2

        if(greater_prob)break; //if a duplicate track is found in a collection with greater Chi^2 probability for Vertex fit then break out of comparison loop
      }//end loop over collection 2 checking
    if(greater_prob)break;//if a duplicate track is found in a collection with greater Chi^2 probability for Vertex fit then one does not need to check the other track the collection will not be stored
  } //end loop over tracks in collection 1
      if(!greater_prob) {
        reco::ConversionRef convRef(convCollHandle, icoll1);
        cleanedConversions->push_back(convRef);
      }
    }//end loop over collection 1

  iEvent.put(cleanedConversions);
}


// ------------ method called once each job just before starting event loop  ------------
void 
ConversionCleanerForPF::beginRun(const edm::Run& run,
             const EventSetup& iSetup) {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ConversionCleanerForPF::endRun(const edm::Run& run,
           const EventSetup& iSetup) {
}
