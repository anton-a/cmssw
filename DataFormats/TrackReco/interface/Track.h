#ifndef TrackReco_Track_h
#define TrackReco_Track_h
/** \class reco::Track Track.h DataFormats/TrackReco/interface/Track.h
 *
 * This class describes the reconstructed tracks that are stored in the AOD and
 * RECO. It also contains a reference to more detailed information about each 
 * track, that is stoed in the TrackExtra object, available only in RECO.
 * 
 * Note that most of the functions provided in this Track class rely on the existance of
 * the TrackExtra object, so will not work on AOD.
 *
 * The most useful functions are those provided in the TrackBase class from which this 
 * inherits, all of which work on AOD.  
 *
 * \author Luca Lista, INFN
 *
 * \version $Id: Track.h,v 1.47.18.1 2013/05/01 15:11:26 anastass Exp $
 *
 */
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h" 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

namespace reco {

  class Track : public TrackBase {
  public:
    /// default constructor
    Track() { }
    /// virtual destructor 
    virtual ~Track();
    /// constructor from fit parameters and error matrix  
    Track( double chi2, double ndof, const Point & referencePoint,
	   const Vector & momentum, int charge, const CovarianceMatrix &,
	  TrackAlgorithm=undefAlgorithm, TrackQuality quality=undefQuality);
    /// constructor from fit parameters, error matrix, and outer position and momentum
    Track( double chi2, double ndof, const Point & referencePoint,
     const Vector & momentum, int charge, const CovarianceMatrix & cov,
     math::XYZPoint innerPosition, math::XYZVector innerMomentum,
     math::XYZPoint outerPosition, math::XYZVector outerMomentum,
     TrackAlgorithm=undefAlgorithm, TrackQuality quality=undefQuality);
    
    /// return true if the outermost hit is valid
    bool outerOk() const { return extra_->outerOk(); }
    /// return true if the innermost hit is valid
    bool innerOk() const { return extra_->innerOk(); }
    
/*    
    /// position of the innermost hit
    const math::XYZPoint & innerPosition()  const { return extra_->innerPosition(); }
    /// momentum vector at the innermost hit position
    const math::XYZVector & innerMomentum() const { return extra_->innerMomentum(); }
*/
    
/*	 
	 /// position of the outermost hit
    const math::XYZPoint & outerPosition()  const { return extra_->outerPosition(); }
    /// momentum vector at the outermost hit position
    const math::XYZVector & outerMomentum() const { return extra_->outerMomentum(); }
*/	 


    /// innermost hit position
    const math::XYZPoint & innerPosition() const { return innerPosition_; }
    /// momentum vector at innermost hit position
    const math::XYZVector & innerMomentum() const { return innerMomentum_; }

	  /// outermost hit position
    const math::XYZPoint & outerPosition() const { return outerPosition_; }
    /// momentum vector at outermost hit position
    const math::XYZVector & outerMomentum() const { return outerMomentum_; }


    /// outermost trajectory state curvilinear errors
    CovarianceMatrix outerStateCovariance() const { return extra_->outerStateCovariance(); }
    /// innermost trajectory state curvilinear errors
    CovarianceMatrix innerStateCovariance() const { return extra_->innerStateCovariance(); }
    /// fill outermost trajectory state curvilinear errors
    CovarianceMatrix & fillOuter( CovarianceMatrix & v ) const { return extra_->fillOuter( v ); }
    /// fill outermost trajectory state curvilinear errors
    CovarianceMatrix & fillInner( CovarianceMatrix & v ) const { return extra_->fillInner( v ); }
    /// DetId of the detector on which surface the outermost state is located
    unsigned int outerDetId() const { return extra_->outerDetId(); }
    /// DetId of the detector on which surface the innermost state is located
    unsigned int innerDetId() const { return extra_->innerDetId(); }
    /// Iterator to first hit on the track.
    trackingRecHit_iterator recHitsBegin() const { return extra_->recHitsBegin(); }
    /// Iterator to last hit on the track.
    trackingRecHit_iterator recHitsEnd() const { return extra_->recHitsEnd(); }
    /// Get i-th hit on the track.
    TrackingRecHitRef recHit( size_t i ) const { return extra_->recHit( i ); }
    /// Get number of RecHits. (Warning, this includes invalid hits, which are not physical hits).
    size_t recHitsSize() const { return extra_->recHitsSize(); }

	 /*
    /// x coordinate of momentum vector at the outermost hit position
    double outerPx()     const { return extra_->outerPx(); }
    /// y coordinate of momentum vector at the outermost hit position
    double outerPy()     const { return extra_->outerPy(); }
    /// z coordinate of momentum vector at the outermost hit position
    double outerPz()     const { return extra_->outerPz(); }
    /// x coordinate of the outermost hit position
    double outerX()      const { return extra_->outerX(); }
    /// y coordinate of the outermost hit position
    double outerY()      const { return extra_->outerY(); }
    /// z coordinate of the outermost hit position
    double outerZ()      const { return extra_->outerZ(); }
    /// magnitude of momentum vector at the outermost hit position
    double outerP()      const { return extra_->outerP(); }
    /// transverse momentum at the outermost hit position
    double outerPt()     const { return extra_->outerPt(); }
    /// azimuthal angle of the  momentum vector at the outermost hit position
    double outerPhi()    const { return extra_->outerPhi(); }
    /// pseudorapidity of the  momentum vector at the outermost hit position
    double outerEta()    const { return extra_->outerEta(); }
    /// polar angle of the  momentum vector at the outermost hit position
    double outerTheta()  const { return extra_->outerTheta(); }    
    /// polar radius of the outermost hit position
    double outerRadius() const { return extra_->outerRadius(); }
	*/


	 // Access information on outer position/momentum directly since they are part of track now (AA)
	 /// x coordinate of momentum vector at the outermost hit position
    double outerPx() const { return outerMomentum_.X(); }
    /// y coordinate of momentum vector at the outermost hit position
    double outerPy() const { return outerMomentum_.Y(); }
    /// z coordinate of momentum vector at the outermost hit position
    double outerPz() const { return outerMomentum_.Z(); }
    /// x coordinate the outermost hit position
    double outerX() const { return outerPosition_.X(); }
    /// y coordinate the outermost hit position
    double outerY() const { return outerPosition_.Y(); }
    /// z coordinate the outermost hit position
    double outerZ() const { return outerPosition_.Z(); }
    /// magnitude of momentum vector at the outermost hit position    
    double outerP() const { return outerMomentum().R(); }
    /// transverse momentum at the outermost hit position    
    double outerPt() const { return outerMomentum().Rho(); }
    /// azimuthal angle of the  momentum vector at the outermost hit position
    double outerPhi() const { return outerMomentum().Phi(); }
    /// pseudorapidity the  momentum vector at the outermost hit position
    double outerEta() const { return outerMomentum().Eta(); }
    /// polar angle of the  momentum vector at the outermost hit position
    double outerTheta() const { return outerMomentum().Theta(); }
    /// polar radius of the outermost hit position
    double outerRadius() const { return outerPosition().Rho(); }

	 
	 /// set the inner position (AA)
	 void setInnerPosition( const math::XYZPoint innerPosition ) {innerPosition_ = innerPosition;}
	 /// set the outer momentum (AA)
	 void setInnerMomentum( const math::XYZVector innerMomentum ) {innerMomentum_ = innerMomentum;}
   /// set the outer position (AA)
   void setOuterPosition( const math::XYZPoint outerPosition ) {outerPosition_ = outerPosition;}
   /// set the outer momentum (AA)
   void setOuterMomentum( const math::XYZVector outerMomentum ) {outerMomentum_ = outerMomentum;}


	 
	  /// set reference to "extra" object
    void setExtra( const TrackExtraRef & ref ) { extra_ = ref; }
    /// reference to "extra" object
    const TrackExtraRef & extra() const { return extra_; }

    /// Number of valid hits on track.
    unsigned short found() const { return  numberOfValidHits(); }
    /// Number of lost (=invalid) hits on track.
    unsigned short lost() const {return  numberOfLostHits();  }

    /// direction of how the hits were sorted in the original seed
    PropagationDirection seedDirection() const {return extra_->seedDirection();}

    /**  return the edm::reference to the trajectory seed in the original
     *   seeds collection. If the collection has been dropped from the
     *   Event, the reference may be invalid. Its validity should be tested,
     *   before the reference is actually used. 
     */
    edm::RefToBase<TrajectorySeed> seedRef() const { return extra_->seedRef(); }

    ///  Access the lightweight track residuals; these are stored in
    ///  TrackExtra and provide residual information with 4 bits of
    ///  precision per hit
    const TrackResiduals &residuals () const { return extra_->residuals(); }
    /// return the residual (local x/y) for the hit in the ith position;
    /// this position is aligned with the position in the HitPattern
    double residualX (int position) const;
    double residualY (int position) const;

  private:
    /// Reference to additional information stored only on RECO.
    TrackExtraRef extra_;

	  /// information (previously in extra) used PF re-reco from AOD (Anton Anastassov)
	 
	  /// innermost hit position
    math::XYZPoint innerPosition_;
    /// momentum vector at innermost hit position
    math::XYZVector innerMomentum_;

    /// outermost hit position
    math::XYZPoint outerPosition_;
    /// momentum vector at outermost hit position
    math::XYZVector outerMomentum_;


  };

}

#endif
