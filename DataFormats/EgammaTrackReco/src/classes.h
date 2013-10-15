#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/CLHEP/interface/Migration.h" 
#include <boost/cstdint.hpp> 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h" 
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateSeedAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/EgammaTrackReco/interface/TrackCandidateSuperClusterAssociation.h"
#include "DataFormats/EgammaTrackReco/interface/TrackSuperClusterAssociation.h"
#include "DataFormats/EgammaTrackReco/interface/TrackCandidateCaloClusterAssociation.h"
#include "DataFormats/EgammaTrackReco/interface/TrackCaloClusterAssociation.h"
#include "DataFormats/Common/interface/RefProd.h" 
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" 
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h" 
#include "DataFormats/EgammaTrackReco/interface/ConversionTrack.h" 
#include "DataFormats/EgammaTrackReco/interface/ConversionTrackFwd.h" 
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EgammaTrackReco/interface/GsfTrackToBremConvTracksMapFwd.h"

namespace {
  struct dictionary {

    reco::TrackCandidateSuperClusterAssociationCollection v5;
    edm::Wrapper<edm::ValueMap<reco::SuperClusterRef> > m5;

    reco::TrackSuperClusterAssociationCollection v6;
    edm::Wrapper<edm::ValueMap<reco::SuperClusterRef> > m6;

    reco::TrackCandidateCaloClusterPtrVectorAssociation v7;
    edm::Wrapper<edm::ValueMap<reco::CaloClusterPtrVector> > m7;

    reco::TrackCandidateCaloClusterPtrAssociation v8;
    edm::Wrapper<edm::ValueMap<reco::CaloClusterPtr> > m8;

    reco::TrackCaloClusterPtrVectorAssociation v9;
    edm::Wrapper<edm::ValueMap<reco::CaloClusterPtrVector> > m9;

    reco::TrackCaloClusterPtrAssociation v10;
    edm::Wrapper<edm::ValueMap<reco::CaloClusterPtr> > m10;

    reco::ConversionTrackCollection v11;
    edm::Wrapper<reco::ConversionTrackCollection> m11;    

    edm::helpers::KeyVal<edm::RefProd<std::vector<reco::GsfTrack> >,edm::RefProd<std::vector<reco::Track> > > am0;
    edm::helpers::KeyVal<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,std::vector<std::pair<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,float> > > amf1;
    edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::GsfTrack>,std::vector<reco::Track>,float,unsigned int> > amf2;
    edm::Wrapper<edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::GsfTrack>,std::vector<reco::Track>,float,unsigned int> > > amf3;
    std::map<unsigned int,edm::helpers::KeyVal<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,std::vector<std::pair<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,float> > > > amf4;
    edm::helpers::KeyVal<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,std::vector<std::pair<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,int> > > am1;
    edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::GsfTrack>,std::vector<reco::Track>,int,unsigned int> > am2;
    edm::Wrapper<edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::GsfTrack>,std::vector<reco::Track>,int,unsigned int> > > am3;
    std::map<unsigned int,edm::helpers::KeyVal<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,std::vector<std::pair<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,int> > > > am4;

    edm::helpers::KeyVal<edm::RefProd<std::vector<reco::Track> >,edm::RefProd<std::vector<reco::GsfTrack> > > ma0;
    edm::helpers::KeyVal<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,std::vector<std::pair<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,int> > > ma1;
    edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::Track>,std::vector<reco::GsfTrack>,int,unsigned int> > ma2;

    edm::Wrapper<edm::AssociationMap<edm::OneToManyWithQuality<std::vector<reco::Track>,std::vector<reco::GsfTrack>,int,unsigned int> > > ma3;
    std::map<unsigned int,edm::helpers::KeyVal<edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> >,std::vector<std::pair<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,int> > > > ma4;

    std::vector<std::pair<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,int> > ma5;
    std::pair<edm::Ref<std::vector<reco::GsfTrack>,reco::GsfTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::GsfTrack>,reco::GsfTrack> >,int> ma6;



  };
}

