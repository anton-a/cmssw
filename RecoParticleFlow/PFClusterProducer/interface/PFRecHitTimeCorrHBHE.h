#ifndef RecoParticleFlow_PFClusterProducer_PFRecHitTimeCorrHBHE_h_
#define RecoParticleFlow_PFClusterProducer_PFRecHitTimeCorrHBHE_h_

#include <cmath>

// Move the corrections once it is decided if they are applied to normal rechits
// and what informaion will be stored in the PFRecHit
// The corrections can be expanded to include time from pulse shape
//
// Anton

/// Linearize the HBHE TDC time on depth-by-depth basis from the measured time
/// Temporary implementation for testing of PF modifications for Phase 2 upgrades.


class PFRecHitTimeCorr {

  public:

  /// return corrected TDC time for HB 
  static float tdcCorrHB(float rawTime, float rhE, int detDepth );
  
  /// return corrected TDC time for HE   
  static float tdcCorrHE(float rawTime, float rhE, int detDepth );

// other corretions
};

#endif



