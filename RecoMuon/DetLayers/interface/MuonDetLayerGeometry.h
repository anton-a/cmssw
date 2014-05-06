#ifndef DetLayers_MuonDetLayerGeometry_h
#define DetLayers_MuonDetLayerGeometry_h

/** \class MuonDetLayerGeometry
 *
 *  Provide access to the DetLayers of muon detectors.
 *
 *  \author N. Amapane - CERN
 */

#include "DataFormats/DetId/interface/DetId.h"
#include "TrackingTools/DetLayers/interface/DetLayerGeometry.h"
#include <vector>
#include <map>

class DetLayer;

class MuonDetLayerGeometry : public DetLayerGeometry{
 public:

  /// Constructor
  MuonDetLayerGeometry();

  friend class MuonDetLayerGeometryESProducer;  

  /// Destructor
  virtual ~MuonDetLayerGeometry();

  /// return the DT DetLayers (barrel), inside-out
  const std::vector<const DetLayer*>& allDTLayers() const;

  /// return the CSC DetLayers (endcap), -Z to +Z
  const std::vector<const DetLayer*>& allCSCLayers() const;

  /// return the forward (+Z) CSC DetLayers, inside-out
  const std::vector<const DetLayer*>& forwardCSCLayers() const;

  /// return the backward (-Z) CSC DetLayers, inside-out
  const std::vector<const DetLayer*>& backwardCSCLayers() const;

  /// return all RPC DetLayers, order: backward, barrel, forward
  const std::vector<const DetLayer*>& allRPCLayers() const;

  /// return the barrel RPC DetLayers, inside-out
  const std::vector<const DetLayer*>& barrelRPCLayers() const;

  /// return the endcap RPC DetLayers, -Z to +Z
  const std::vector<const DetLayer*>& endcapRPCLayers() const;

  /// return the forward (+Z) RPC DetLayers, inside-out
  const std::vector<const DetLayer*>& forwardRPCLayers() const;

  /// return the backward (-Z) RPC DetLayers, inside-out
  const std::vector<const DetLayer*>& backwardRPCLayers() const;

  /// return all layers (DT+CSC+RPC), order: backward, barrel, forward
  const std::vector<const DetLayer*>& allLayers() const;

  /// return all barrel DetLayers (DT+RPC), inside-out
  const std::vector<const DetLayer*>& allBarrelLayers() const;

  /// return all endcap DetLayers (CSC+RPC), -Z to +Z
  const std::vector<const DetLayer*>& allEndcapLayers() const;

  /// return all forward (+Z) layers (CSC+RPC), inside-out
  const std::vector<const DetLayer*>& allForwardLayers() const;

  /// return all backward (-Z) layers (CSC+RPC), inside-out
  const std::vector<const DetLayer*>& allBackwardLayers() const;
  
  /// return the DetLayer which correspond to a certain DetId
  virtual const DetLayer* idToLayer(const DetId& detId) const override;

 private:
  /// Add CSC layers 
  /// csclayers.first=forward (+Z), csclayers.second=backward (-Z)
  /// both vectors are ASSUMED to be sorted inside-out
  void addCSCLayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& csclayers);

  //. Add DT layers; dtlayers is ASSUMED to be sorted inside-out
  void addDTLayers(const std::vector<DetLayer*>& dtlayers);

  /// Add RPC layers
  /// endcapRPCLayers.first=forward (+Z), endcapRPCLayers.second=backward (-Z)
  /// All three vectors are ASSUMED to be sorted inside-out
  void addRPCLayers(const std::vector<DetLayer*>& barrelRPCLayers, const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& endcapRPCLayers);

  
  DetId makeDetLayerId(const DetLayer* detLayer) const;
  
  void sortLayers();

  std::vector<const DetLayer*> cscLayers_fw;
  std::vector<const DetLayer*> cscLayers_bk;
  std::vector<const DetLayer*> cscLayers_all;
  std::vector<const DetLayer*> rpcLayers_all;
  std::vector<const DetLayer*> rpcLayers_endcap;
  std::vector<const DetLayer*> rpcLayers_fw;
  std::vector<const DetLayer*> rpcLayers_bk;
  std::vector<const DetLayer*> rpcLayers_barrel;
  std::vector<const DetLayer*> dtLayers;
  std::vector<const DetLayer*> allForward;
  std::vector<const DetLayer*> allBackward;
  std::vector<const DetLayer*> allEndcap;
  std::vector<const DetLayer*> allBarrel;
  std::vector<const DetLayer*> allDetLayers;
    
  std::map<DetId,const DetLayer*> detLayersMap;
};
#endif

