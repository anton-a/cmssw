#include <memory>
#include "RecoParticleFlow/PFTracking/interface/ConvBremGsfTrackMapProducer.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoParticleFlow/PFTracking/interface/ConvBremTrackFinder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// the new map of gsftracks to brem conv tracks stored in the event (AA)
#include "DataFormats/EgammaTrackReco/interface/GsfTrackToBremConvTracksMapFwd.h"

using namespace std;
using namespace edm;
using namespace reco;
ConvBremGsfTrackMapProducer::ConvBremGsfTrackMapProducer(const ParameterSet& iConfig):
  pfTransformer_(0),
  convBremTrackFinder_(0)
{
  produces<reco::GsfTrackToBremConvTracksMap>();
  
  trackCollection_ = consumes<TrackCollection>(iConfig.getParameter<InputTag>
    ("TrackCollection"));
  
  useQuality_ = iConfig.getParameter<bool>
    ("UseQuality");
  
  gsfTrackLabel_ = consumes<GsfTrackCollection>(iConfig.getParameter<InputTag>
    ("GsfTrackModuleLabel"));  

  pfEcalClusters_ = consumes<PFClusterCollection>(iConfig.getParameter<InputTag>
    ("PFEcalClusters"));

  trackQuality_ = reco::TrackBase::qualityByName(iConfig.getParameter<std::string>
						 ("TrackQuality"));

  vtx_h = consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>
    ("PrimaryVertexLabel"));

  // set parameter for convBremFinder
  useConvBremFinder_ = iConfig.getParameter<bool>
    ("useConvBremFinder");

  mvaConvBremFinderID_ = iConfig.getParameter<double>
    ("pf_convBremFinderID_mvaCut");
  
  string mvaWeightFileConvBrem = iConfig.getParameter<string>
    ("pf_convBremFinderID_mvaWeightFile");
  
  
  if(useConvBremFinder_) 
    path_mvaWeightFileConvBrem_ = edm::FileInPath ( mvaWeightFileConvBrem.c_str() ).fullPath();


}

ConvBremGsfTrackMapProducer::~ConvBremGsfTrackMapProducer()
{
  delete pfTransformer_;
  delete convBremTrackFinder_;
}

void
ConvBremGsfTrackMapProducer::produce(Event& iEvent, const EventSetup& iSetup)
{
  
  bool debug = false;

  // Crate a collection to hold the GSF track - conv brem track association
   auto_ptr< reco::GsfTrackToBremConvTracksMap > 
   gsfTrackToConvTracksMap(new reco::GsfTrackToBremConvTracksMap);
 
   Handle<reco::TrackCollection> theTrackCollection;
   iEvent.getByToken(trackCollection_, theTrackCollection);
 
  //read track collection
  Handle<GsfTrackCollection> gsftrackcoll;
  iEvent.getByToken(gsfTrackLabel_,gsftrackcoll);
  GsfTrackCollection gsftracks = *(gsftrackcoll.product());

  //PFClusters
  Handle<PFClusterCollection> theECPfClustCollection;
  iEvent.getByToken(pfEcalClusters_,theECPfClustCollection);
  const PFClusterCollection& theEcalClusters = *(theECPfClustCollection.product());

  //Primary Vertexes
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByToken(vtx_h,thePrimaryVertexColl);

  if (debug) cout << " ######## ConvBremGsfTrackMaptProdcuer:Entering  " << endl;

  for (unsigned int igsf=0; igsf<gsftracks.size();igsf++) {
    
    if (debug)
      cout << " NEW GSF: pt " << gsftracks[igsf].ptMode() 
	    << " eta " << gsftracks[igsf].etaMode() << " phi " <<   gsftracks[igsf].phiMode() << endl; 

    if (convBremTrackFinder_->foundConvBremTrack(theTrackCollection, thePrimaryVertexColl,
                                                 theEcalClusters, gsftracks[igsf])) {

       const vector<TrackRef>& convBremTracks(convBremTrackFinder_->getConvBremTracks());
       const vector<float>& mvaScores =  convBremTrackFinder_->getConvBremTrackMVAScores();
       
       if(debug) cout << " ######## ConvBremGsfTrackMaptProdcuer:foundConvBremTrack  " << endl;

       // Note that we make an entry in the map only when conv brem tracks were found
       // When acessing the map one should check if there is a an entry for the particular gsftrack (AA)
       for (unsigned int itk = 0; itk < convBremTracks.size(); ++itk) {
       
          if (debug) cout << " ConvBremTrack Pt: " << convBremTracks[itk]->pt() 
	            << "\teta, phi: " << convBremTracks[itk]->eta() << ", " << convBremTracks[itk]->phi() 
              << "\tMVA score: " << mvaScores[itk] << endl;

          // without quality of match
//        gsfTrackToConvTracksMap->insert(reco::GsfTrackRef(gsftrackcoll, igsf), convBremTracks[itk]);

          // with quality of matching
          gsfTrackToConvTracksMap->insert(reco::GsfTrackRef(gsftrackcoll, igsf), 
                                        std::make_pair(convBremTracks[itk], mvaScores[itk]) );


     } // loop over converted brem tracks
  
    } // if found conv brem tracks

  } // end loop on the GSF track
  
  
  if (debug) cout << "\n Before putting gsfTrackToConvTracksMap in the event\n";
 // out the map in the event
 iEvent.put(gsfTrackToConvTracksMap);

}

// ------------ method called once each job just before starting event loop  ------------
void 
ConvBremGsfTrackMapProducer::beginRun(const edm::Run& run,
			  const EventSetup& iSetup)
{
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);


  ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());


  pfTransformer_= new PFTrackTransformer(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));
  pfTransformer_->OnlyProp();

  mtsTransform_ = MultiTrajectoryStateTransform(tracker.product(),magneticField.product());

  convBremTrackFinder_ = new ConvBremTrackFinder(thebuilder,
						 mvaConvBremFinderID_,
						 path_mvaWeightFileConvBrem_,
						 *pfTransformer_,
						 mtsTransform_,
						 trackQuality_);

}  // beginRun

// ------------ method called once each job just after ending the event loop  ------------
void 
ConvBremGsfTrackMapProducer::endRun(const edm::Run&,const edm::EventSetup&) {
  delete pfTransformer_;
  delete convBremTrackFinder_;
  pfTransformer_ = nullptr;
  convBremTrackFinder_ = nullptr;
}
