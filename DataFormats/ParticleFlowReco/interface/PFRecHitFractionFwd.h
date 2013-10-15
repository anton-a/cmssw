#ifndef ParticleFlowReco_PFRecHitFractionFwd_h
#define ParticleFlowReco_PFRecHitFractionFwd_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace reco {
  class PFRecHitFraction;

  /// collection of PFRecHitFraction objects
  typedef std::vector<std::vector<PFRecHitFraction> > PFRecHitFractionCollection;

  /// persistent reference to PFRecHitFraction objects
  typedef edm::Ref<PFRecHitFractionCollection> PFRecHitFractionRef;

  /// reference to PFRecHitFraction collection
  typedef edm::RefProd<PFRecHitFractionCollection> PFRecHitFractionRefProd;

  /// vector of references to PFRecHitFraction objects all in the same collection
  typedef edm::RefVector<PFRecHitFractionCollection> PFRecHitFractionRefVector;

  /// iterator over a vector of references to PFRecHitFraction objects
  typedef PFRecHitFractionRefVector::iterator PFRecHitFraction_iterator;
}

#endif
