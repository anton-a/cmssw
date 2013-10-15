#ifndef DataFormats_EgammaTrackReco_GsfTrackToConvTracksMapFwd_h
#define DataFormats_EgammaTrackReco_GsfTrackToConvTracksMapFwd_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/Common/interface/AssociationMap.h"

/// Container for association of each GSF track to brem conversion tracks

namespace reco {


// with quality of brem conv track association (AA)
   typedef edm::AssociationMap<edm::OneToManyWithQuality<reco::GsfTrackCollection, reco::TrackCollection, float> > GsfTrackToBremConvTracksMap;


/*
   // previous version: without quality
   typedef edm::AssociationMap<edm::OneToMany<reco::GsfTrackCollection, reco::TrackCollection> > GsfTrackToBremConvTracksMap;
   typedef  GsfTrackToBremConvTracksMap::value_type GsfTrackToBremConvTracks;
 
   /// reference to an object in a collection of GsfTrackToBremConvTracksMap objects
   typedef edm::Ref<GsfTrackToBremConvTracksMap> GsfTrackToBremConvTracksMapRef;

   /// reference to a collection of GsfTrackToBremConvTracksMap object
   typedef edm::RefProd<GsfTrackToBremConvTracksMap> GsfTrackToBremConvTracksMapRefProd;

   /// vector of references to objects in the same colletion of GsfTrackToBremConvTracksMap objects
   typedef edm::RefVector<GsfTrackToBremConvTracksMap> GsfTrackToBremConvTracksMapRefVector;
*/

}

#endif
