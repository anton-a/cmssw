#ifndef ConversionCleanerForPF_H
#define ConversionCleanerForPF_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

/**\class ConversionCleanerForPF 
\brief Cleans the allConversions collection for input to the PFConversionsProducer 

This producer makes a collections of references to the conversions in the 
allConversions collection that pass the requirements for PF.

\author Anton Anastassov
\date   August 2013
*/

class ConversionCleanerForPF : public edm::EDProducer {
public:

  ///Constructor
  explicit ConversionCleanerForPF(const edm::ParameterSet&);
  
  ///Destructor
  ~ConversionCleanerForPF();

private:

  virtual void beginRun(const edm::Run&,const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&,const edm::EventSetup&) override;

  ///Produce collection of references to accepted conversions
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

  ///Input conversion collection
  edm::EDGetTokenT<reco::ConversionCollection> inputConversions_;

};
#endif
