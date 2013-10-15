#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <iostream>

#include "TMVA/Reader.h"

class PFEnergyCalibration;

class ConvBremTrackFinder {
  
 public:
  ConvBremTrackFinder(const TransientTrackBuilder& builder,
		      double mvaBremConvCut,
		      std::string mvaWeightFileConvBrem,
		      const PFTrackTransformer& pfTkTransformer,
		      const MultiTrajectoryStateTransform& mtjstate,
		      const reco::TrackBase::TrackQuality& trackQuality);
  ~ConvBremTrackFinder();
  
  bool foundConvBremTrack(const edm::Handle<reco::TrackCollection>& theTrackColl,
			  const edm::Handle<reco::VertexCollection>& primaryVertex,
			  const reco::PFClusterCollection & theEClus,
			  reco::GsfTrack gsftrack)
  {
    found_ = false;
    runConvBremFinder(theTrackColl,primaryVertex,
		      theEClus,gsftrack);
    return found_;};
  
  
  const std::vector<reco::TrackRef>& getConvBremTracks() const {return  trackRef_vec_;}
  // the same order as the tracks
  const std::vector<float> & getConvBremTrackMVAScores() const {return mvaScores_;}

 private:
  void runConvBremFinder(const edm::Handle<reco::TrackCollection>& theTrackColl,
			 const edm::Handle<reco::VertexCollection>& primaryVertex,
			 const reco::PFClusterCollection & theEClus,
			 reco::GsfTrack gsftrack);


  bool found_;
  TransientTrackBuilder builder_;
  double mvaBremConvCut_;
  std::string mvaWeightFileConvBrem_;
  TMVA::Reader    *tmvaReader_;

  // TMVA input variables
  float secR,secPout,ptRatioGsfKF,sTIP,Epout,detaBremKF,secPin;
  float nHITS1;

  // vector of matched conversion tracks
  std::vector<reco::TrackRef> trackRef_vec_;
  // corresponding (same ordering) of MVA scores
  std::vector<float> mvaScores_;
  
  PFTrackTransformer pfTkTransformer_;     
  PFEnergyCalibration* pfcalib_;
  MultiTrajectoryStateTransform mtjstate_;
  reco::TrackBase::TrackQuality trackQuality_;
};
