#ifndef ConvBremGsfTrackMapProducer_H
#define ConvBremGsfTrackMapProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

/// \brief Abstract
/*!
\author Daniele Benedetti
\author Anton Anastassov
\date April 2013
  ConvBremGsfTrackMapProducer produces a map between KF tracks from converted brems and a GSF track. 
*/


class PFTrackTransformer;
class ConvBremTrackFinder;

class ConvBremGsfTrackMapProducer : public edm::EDProducer {
public:
  
  ///Constructor
  explicit ConvBremGsfTrackMapProducer(const edm::ParameterSet&);
  
  ///Destructor
  ~ConvBremGsfTrackMapProducer();
  
private:
  virtual void beginRun(const edm::Run&,const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&,const edm::EventSetup&) override;
  
  ///Produce the Track collection
  virtual void produce(edm::Event&, const edm::EventSetup&);
  

  /// Inputs collections
  edm::EDGetTokenT<reco::TrackCollection> trackCollection_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTrackLabel_;  
  edm::EDGetTokenT<reco::VertexCollection> vtx_h;
  edm::EDGetTokenT<reco::PFClusterCollection> pfEcalClusters_;


  /// track quality
  bool useQuality_;
  reco::TrackBase::TrackQuality trackQuality_;

  /// Conv Brem Finder
  bool useConvBremFinder_;
  double mvaConvBremFinderID_;
  std::string path_mvaWeightFileConvBrem_;

  /// Transformers  
  PFTrackTransformer *pfTransformer_; 
  MultiTrajectoryStateTransform  mtsTransform_;
  ConvBremTrackFinder *convBremTrackFinder_;

};
#endif
