#ifndef DataFormats_ParticleFlowReco_PFCluster_h
#define DataFormats_ParticleFlowReco_PFCluster_h

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "Math/GenVector/PositionVector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Rtypes.h" 

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFractionFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include <iostream>
#include <vector>



class PFClusterAlgo;

namespace reco {

  /**\class PFCluster
     \brief Particle flow cluster, see clustering algorithm in PFClusterAlgo

     A particle flow cluster is defined by its energy and position, which are 
     calculated from a vector of PFRecHitFraction. This calculation is 
     performed in PFClusterAlgo.

     \todo Clean up this class to a common base (talk to Paolo Meridiani)
     the extra internal stuff (like the vector of PFRecHitFraction's)
     could be moved to a PFClusterExtra.

     \todo Now that PFRecHitFraction's hold a reference to the PFRecHit's, 
     put back the calculation of energy and position to PFCluster. 


     \todo Add an operator+=

     \author Colin Bernet
     \date   July 2006
  */
  class PFCluster : public CaloCluster {
  public:


    typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;
  
    PFCluster() : CaloCluster(CaloCluster::particleFlow), color_(1) {}

    /// constructor
    PFCluster(PFLayer::Layer layer, double energy,
	      double x, double y, double z );



    /// resets clusters parameters
    void reset();
    
    /// add a given fraction of the rechit
    void addRecHitFraction( const reco::PFRecHitFraction& frac);
    
    /// Get the transient PFRecHitFractions. Needed for storing them in a separate collection. Used only in the cluster producer. 
    std::vector< reco::PFRecHitFraction > transientPFRecHitFractions() { return rechits_;}

    /// vector of rechit fractions
    // in the previous format 
    // this accessor is used by consumers of the info, but not internally
    // convert it to use the new container (AA)
//    const std::vector< reco::PFRecHitFraction >& recHitFractions() const 
//      { return rechits_; }

   // use the ref to the separately stored collection (AA) 
    const PFRecHitFractionRef & recHitFractions() const { return recHitFracs_; }

    /// Set the reference to the PFRecHitFractions (AA)
    void setPFRecHitFractionsRef(PFRecHitFractionRef rhf) { recHitFracs_ = rhf; }



    /// set layer
    void setLayer( PFLayer::Layer layer);
    
    /// cluster layer, see PFLayer.h in this directory
    PFLayer::Layer  layer() const;     
    
    /// cluster energy
    double        energy() const {return energy_;}
    
    /// cluster position: rho, eta, phi
    const REPPoint&       positionREP() const {return posrep_;}
    
    /// computes posrep_ once and for all
    void calculatePositionREP() {
      posrep_.SetCoordinates( position_.Rho(), 
			      position_.Eta(), 
			      position_.Phi() ); 
    }
    
    /// \todo move to PFClusterTools
    static double getDepthCorrection(double energy, bool isBelowPS = false,
				     bool isHadron = false);
    
    /// set cluster color (for the PFRootEventManager display)
    void         setColor(int color) {color_ = color;}
    
    /// \return color
    int          color() const {return color_;}
    
    
    PFCluster& operator=(const PFCluster&);
    
    friend    std::ostream& operator<<(std::ostream& out, 
				       const PFCluster& cluster);
    /// counter
    static unsigned     instanceCounter_;
    
    /// \todo move to PFClusterTools
    static void setDepthCorParameters(int mode, 
				      double a, double b, 
				      double ap, double bp ) {
      depthCorMode_ = mode;
      depthCorA_ = a; 
      depthCorB_ = b; 
      depthCorAp_ = ap; 
      depthCorBp_ = bp; 
    } 
    

    /// some classes to make this fit into a template footprint
    /// for RecoPFClusterRefCandidate so we can make jets and MET
    /// out of PFClusters.
    
    /// dummy charge
    double charge() const { return 0;}

    /// transverse momentum, massless approximation
    double pt() const { 
      return (energy() * sin(position_.theta()));
    }

    /// angle
    double theta() const { 
      return position_.theta();
    }
    
    /// dummy vertex access
    math::XYZPoint const & vertex() const { 
      static math::XYZPoint dummyVtx(0,0,0);
      return dummyVtx;      
    }
    double vx() const { return vertex().x(); }
    double vy() const { return vertex().y(); }
    double vz() const { return vertex().z(); }    

  private:
    
    /// vector of rechit fractions (transient)
    // use these for internal manipulation (AA)
    std::vector< reco::PFRecHitFraction >  rechits_;
    
    /// link to a vector of rechit fractions (AA)
    // these will be persistent in RECO 
    PFRecHitFractionRef  recHitFracs_;

    /// cluster position: rho, eta, phi (transient)
    REPPoint            posrep_;
    
    /// \todo move to PFClusterTools
    static int    depthCorMode_;
    
    /// \todo move to PFClusterTools
    static double depthCorA_;
    
    /// \todo move to PFClusterTools
    static double depthCorB_ ;
    
    /// \todo move to PFClusterTools
    static double depthCorAp_;
    
    /// \todo move to PFClusterTools
    static double depthCorBp_;
    
    
    /// color (transient)
    int                 color_;
    
    friend class ::PFClusterAlgo;

  };
}

#endif
