#include <memory>
#include "RecoParticleFlow/PFTracking/plugins/PFConversionProducer.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversionFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

typedef  std::multimap<unsigned, std::vector<unsigned> > BlockMap;
using namespace std;
using namespace edm;


PFConversionProducer::PFConversionProducer(const ParameterSet& iConfig):
  pfTransformer_(0)
{
  produces<reco::PFRecTrackCollection>();
  produces<reco::PFConversionCollection>();

  inputConversionsTag_ = 
    iConfig.getParameter< InputTag >("inputConversions");
  inputVerticesTag_=iConfig.getParameter<edm::InputTag>("inputPrimaryVertices");
}

PFConversionProducer::~PFConversionProducer()
{
  delete pfTransformer_;
}

void
PFConversionProducer::produce(Event& iEvent, const EventSetup& iSetup)
{
  
  //create the empty collections 
  auto_ptr< reco::PFConversionCollection > 
    pfConversionColl (new reco::PFConversionCollection);
  auto_ptr< reco::PFRecTrackCollection > 
    pfRecTrackColl (new reco::PFRecTrackCollection);
  
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());
  reco::PFRecTrackRefProd pfTrackRefProd = iEvent.getRefBeforePut<reco::PFRecTrackCollection>();

/*
  Handle<reco::ConversionCollection> convCollH;
  iEvent.getByLabel(inputConversionsTag_, convCollH);
  const reco::ConversionCollection& convColl = *(convCollH.product());
*/

  Handle<reco::ConversionRefVector> convCollH;
  iEvent.getByLabel(inputConversionsTag_, convCollH);
  const reco::ConversionRefVector& convColl = *(convCollH.product());

  
  Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(inputVerticesTag_, vertexHandle);
  //Find PV for IP calculation, if there is no PV in collection than use dummy 
  reco::Vertex dummy;
  const reco::Vertex* pv=&dummy;    
  if (vertexHandle.isValid()) 
    {
      pv = &*vertexHandle->begin();
    } 
  else 
    { // create a dummy PV
      reco::Vertex::Error e;
      e(0, 0) = 0.0015 * 0.0015;
      e(1, 1) = 0.0015 * 0.0015;
      e(2, 2) = 15. * 15.;
      reco::Vertex::Point p(0, 0, 0);
      dummy = reco::Vertex(p, e, 0, 0, 0);   
    } 
  
  int idx = 0; //index of track in PFRecTrack collection 
  
//multimap<unsigned int, unsigned int> trackmap; //Map of Collections and tracks -> not used? (AA)

  // Fill the collections
  for(unsigned iColl=0; iColl<convColl.size(); iColl++) {
      std::vector<reco::PFRecTrackRef> pfRecTkcoll;	
      
      std::vector<edm::RefToBase<reco::Track> > tracksRefColl = convColl[iColl]->tracks();	  
      // convert the secondary tracks
      for(unsigned it = 0; it < tracksRefColl.size(); it++)
	{
	  reco::TrackRef trackRef = (tracksRefColl[it]).castTo<reco::TrackRef>();      
	  reco::PFRecTrack pfRecTrack( trackRef->charge(), 
				       reco::PFRecTrack::KF, 
				       trackRef.key(), 
				       trackRef );             
	  //std::cout<<"Track Pt "<<trackRef->pt()<<std::endl;
	  Trajectory FakeTraj;
	  bool valid = pfTransformer_->addPoints( pfRecTrack, *trackRef, FakeTraj);
	  if(valid) 
	    {
	      double stip=-999;
	      const reco::PFTrajectoryPoint& atECAL=pfRecTrack.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance);
	      //if extrapolation to ECAL is valid then calculate STIP
	      if(atECAL.isValid())
		{
		  GlobalVector direction(pfRecTrack.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance).position().x(),
					 pfRecTrack.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance).position().y(), 
					 pfRecTrack.extrapolatedPoint(reco::PFTrajectoryPoint::ECALEntrance).position().z());
		  stip = IPTools::signedTransverseImpactParameter(thebuilder.build(*trackRef), direction, *pv).second.significance();
		}
	      pfRecTrack.setSTIP(stip);	    
	      pfRecTkcoll.push_back(reco::PFRecTrackRef( pfTrackRefProd, idx++));    	  
	      pfRecTrackColl->push_back(pfRecTrack);	    
	    }
	}//end loop over tracks
      //store reference to the Conversion collection
      pfConversionColl->push_back( reco::PFConversion(convColl[iColl] , pfRecTkcoll ));
      
      
    }//end loop over collections
  iEvent.put(pfRecTrackColl);
  iEvent.put(pfConversionColl);    
}
  
// ------------ method called once each job just before starting event loop  ------------
void 
PFConversionProducer::beginRun(const edm::Run& run,
			       const EventSetup& iSetup)
{
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  pfTransformer_= new PFTrackTransformer(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));
  pfTransformer_->OnlyProp();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFConversionProducer::endRun(const edm::Run& run,
			     const EventSetup& iSetup) {
  delete pfTransformer_;
  pfTransformer_=nullptr;
}
