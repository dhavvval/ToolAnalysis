#ifndef DummyTriggerHitsOverlay_H
#define DummyTriggerHitsOverlay_H

#include <string>
#include <vector>
#include <map>
#include <cstdint>

#include "Tool.h"
#include "Hit.h"
#include "Geometry.h"

/**
 * \class DummyTriggerHitsOverlay
 *
 * Adds the tank hits of one dummy-trigger (random-trigger) data window to every MC event,
 * after PhaseIIADCHitFinder (BuildWaveformMCHits 1) and before ClusterFinder. See README.md.
 *
 * $Author: D. Ajana $
 * $Date: 2026/10/30 $
 */
class DummyTriggerHitsOverlay: public Tool {

 public:

  DummyTriggerHitsOverlay(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Reads the dummy-trigger windows and the data-to-MC channel map
  bool Execute(); ///< Appends the hits of one window to the MC hit store
  bool Finalise(); ///< Prints the overlay summary

  static constexpr int kOverlayDirectParent = -3; ///< DirectParents code of an overlay hit (-1 dark noise, -2 track not saved)

 private:

  bool LoadTriggers();
  bool BuildChannelMap();
  size_t PickTrigger(const std::string& mcfile, uint64_t mcevent) const;

  int verbosity = 1;
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;

  std::string fOverlayFile;
  std::string fTreeName = "Event";
  std::string fHitStoreName = "WaveformMCHits";
  int fSeed = 0;
  double fTimeOffset = 0.;

  Geometry* geom = nullptr;

  // one entry per dummy trigger with tank data; hits of trigger i are [fTrigFirstHit[i], fTrigFirstHit[i+1])
  std::vector<long> fTrigEntry;
  std::vector<uint64_t> fTrigTankTime;
  std::vector<size_t> fTrigFirstHit;
  std::vector<int> fHitChankey; // data channel key
  std::vector<double> fHitT;
  std::vector<double> fHitPE;

  std::map<unsigned long, unsigned long> fDataToMCChankey;

  long fNEvents = 0;
  long fNAdded = 0;
  long fNUnmapped = 0;
  long fNChannelOff = 0;
  std::vector<int> fTimesUsed;

};


#endif
