#include "RecoParticleFlow/PFClusterProducer/interface/PFRecHitTimeCorrHBHE.h"

float PFRecHitTimeCorr::tdcCorrHB(float rawTime, float rhE, int detDepth ) {

  if (rhE > 40.0) rhE = 40.0;

  float p0 = 0.0;
  float p1 = 0.0;
  float p2 = 0.0;
  float p3 = 0.0;
  float p4 = 0.0;
  
  switch (detDepth) {
    case 1: 
    p0 =         3.28587e+000;
    p1 =        -6.75780e-001;
    p2 =         7.33208e+000;
    p3 =        -1.60925e-001;
    p4 =         1.58950e-003;
    break;
    case 2:
    p0 =         3.18140e+000;
    p1 =        -6.35764e-001;
    p2 =         6.84248e+000;
    p3 =        -1.77508e-001;
    p4 =         1.86179e-003;
    break;
    case 3:
    p0 =         3.25434e+000;
    p1 =        -6.62655e-001;
    p2 =         6.03523e+000;
    p3 =        -2.58276e-001;
    p4 =         3.86700e-003;
    break;
  }
    
  return rawTime - (exp(p0+p1*rhE) + p2 + p3*rhE +p4*rhE*rhE);
    
}

float PFRecHitTimeCorr::tdcCorrHE(float rawTime, float rhE, int detDepth ) {

  if (rhE > 40.0) rhE = 40;

  float p0 = 0.0;
  float p1 = 0.0;
  float p2 = 0.0;
  float p3 = 0.0;
  float p4 = 0.0;
  
  switch (detDepth) {
    case 1:
    p0 =         3.44969e+000;
    p1 =        -7.74797e-001;
    p2 =         7.07454e+000;
    p3 =        -1.75267e-001;
    p4 =         2.08585e-003;
    break;
    
    case 2:
    p0 =         3.42694e+000;
    p1 =        -7.38052e-001;
    p2 =         7.32696e+000;
    p3 =        -2.14835e-001;
    p4 =         2.98468e-003;
    break;
    
    case 3: 
    p0 =         3.44122e+000;
    p1 =        -7.32644e-001;
    p2 =         6.68813e+000;
    p3 =        -2.13227e-001;
    p4 =         2.89048e-003;
    break;
    
    case 4: 
    p0 =         3.44769e+000;
    p1 =        -7.36780e-001;
    p2 =         5.77368e+000;
    p3 =        -2.06017e-001;
    p4 =         2.72580e-003;
    break;
    
    case 5:
    p0 =         3.45880e+000;
    p1 =        -7.49052e-001;
    p2 =         4.39402e+000;
    p3 =        -2.16175e-001;
    p4 =         2.65140e-003;
    break;
  }

  return rawTime - (exp(p0+p1*rhE) + p2 + p3*rhE +p4*rhE*rhE);

}