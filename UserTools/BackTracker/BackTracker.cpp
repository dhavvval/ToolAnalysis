#include "BackTracker.h"
#include "ANNIEconstants.h"

BackTracker::BackTracker():Tool(){}

// To sort
struct sort_by_charge {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};


bool BackTracker::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load my config parameters
  bool gotVerbosity = m_variables.Get("verbosity",verbosity);
  if (!gotVerbosity) {
    verbosity = 0;
    logmessage = "BackTracker::Initialize: \"verbosity\" not set in the config, defaulting to 0";
    Log(logmessage, v_error, verbosity);
  }

  bool gotUsePulseWindowMatching = m_variables.Get("UsePulseWindowMatching", fUsePulseWindowMatching);
  if (!gotUsePulseWindowMatching) fUsePulseWindowMatching = true;


  // Set up the pointers we're going to save. No need to 
  // delete them at Finalize, the store will handle it
  fClusterToBestParticleID  = new std::map<double, int>;
  fClusterToBestParticlePDG = new std::map<double, int>;
  fClusterEfficiency        = new std::map<double, double>;
  fClusterPurity            = new std::map<double, double>;
  fClusterTotalCharge       = new std::map<double, double>;
  fHitToDirectParents       = new std::map<unsigned long, std::map<double, std::vector<int>>>;
  //fClusterHitToDirectParentTrackIDs = new std::map<double, std::vector<std::vector<int>>>;
  
  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Execute()
{
  if (!LoadFromStores())
    return false;

  if (fUsePulseWindowMatching) {
    // Required tool order: PMTWaveformSim -> PhaseIIADCHitFinder -> BackTracker
    DirectParentsFromClockTickWindows();
  }

  fClusterToBestParticleID ->clear();
  fClusterToBestParticlePDG->clear();
  fClusterEfficiency       ->clear();
  fClusterPurity           ->clear();
  fClusterTotalCharge      ->clear();
  fParticleToTankTotalCharge.clear();
  SumParticleTankCharge();

  // Loop over the clusters and do the things
  for (std::pair<double, std::vector<MCHit>>&& apair : *fClusterMapMC) {
    int prtId = -5;
    int prtPdg = -5;
    double eff = -5;
    double pur = -5;
    double totalCharge = 0;

    MatchMCParticle(apair.second, prtId, prtPdg, eff, pur, totalCharge);

    fClusterToBestParticleID ->emplace(apair.first, prtId);
    fClusterToBestParticlePDG->emplace(apair.first, prtPdg);
    fClusterEfficiency       ->emplace(apair.first, eff);
    fClusterPurity           ->emplace(apair.first, pur);
    fClusterTotalCharge      ->emplace(apair.first, totalCharge);

  }

  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticleID",  fClusterToBestParticleID );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  m_data->Stores.at("ANNIEEvent")->Set("ClusterEfficiency",        fClusterEfficiency       );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterPurity",            fClusterPurity           );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterTotalCharge",       fClusterTotalCharge      );
  m_data->Stores.at("ANNIEEvent")->Set("HitToDirectParents",       fHitToDirectParents      );

  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
void BackTracker::SumParticleTankCharge()
{
  for (auto mcHitsIt : *fMCHitsMap) {
    std::vector<MCHit> mcHits = mcHitsIt.second;
    for (uint mcHitIdx = 0; mcHitIdx < mcHits.size(); ++mcHitIdx) {

      // technically a MCHit could have multiple parents, but they don't appear to in practice
      // skip any cases we come across
      std::vector<int> parentIdxs = *(mcHits[mcHitIdx].GetParents());
      if (parentIdxs.size() != 1) continue;
      
      int particleId = -5;
      for (auto it : *fMCParticleIndexMap) {
	if (it.second == parentIdxs[0]) particleId = it.first;
      }
      if (particleId == -5) continue;
	
      double depositedCharge = mcHits[mcHitIdx].GetCharge();      
      if (!fParticleToTankTotalCharge.count(particleId)) 
	fParticleToTankTotalCharge.emplace(particleId, depositedCharge);
      else 
	fParticleToTankTotalCharge.at(particleId) += depositedCharge;
    }    
  }
}

//------------------------------------------------------------------------------
void BackTracker::MatchMCParticle(std::vector<MCHit> const &mchits, int &prtId, int &prtPdg, double &eff, double &pur, double &totalCharge)
{
  // Loop over the hits and get all of their parents and the energy that each one contributed
  //  be sure to bunch up all neutronic contributions
  std::map<int, double> mapParticleToTotalClusterCharge;
  totalCharge = 0;

  for (auto mchit : mchits) {    
    std::vector<int> parentIdxs = *(mchit.GetParents());
    if (parentIdxs.size() != 1) {
      logmessage = "BackTracker::MatchMCParticle: this MCHit has ";
      logmessage += std::to_string(parentIdxs.size()) + " parents!";
      Log(logmessage, v_debug, verbosity);
      continue;
    }
    
    int particleId = -5;
    for (auto it : *fMCParticleIndexMap) {
      if (it.second == parentIdxs[0]) particleId = it.first;
    }
    if (particleId == -5) continue;
    
    double depositedCharge = mchit.GetCharge();
    totalCharge += depositedCharge;
    
    if (mapParticleToTotalClusterCharge.count(particleId) == 0) 
      mapParticleToTotalClusterCharge.emplace(particleId, depositedCharge);
    else
      mapParticleToTotalClusterCharge[particleId] += depositedCharge;    
  }       

  // Loop over the particleIds to find the primary contributer to the cluster
  double maxCharge = 0;
  for (auto apair : mapParticleToTotalClusterCharge) {
    if (apair.second > maxCharge) {
      maxCharge = apair.second;
      prtId = apair.first;
    }
  }

  // Check that we have some charge, if not then something is wrong so pass back all -5
  if (totalCharge > 0) {
    eff = maxCharge/fParticleToTankTotalCharge.at(prtId);
    pur = maxCharge/totalCharge;
    prtPdg = (fMCParticles->at(fMCParticleIndexMap->at(prtId))).GetPdgCode();
  } else {
    prtId = -5;
    eff = -5;
    pur = -5;
    totalCharge = -5;
  }

  logmessage = "BackTracker::MatchMCParticle: best particleId is : ";
  logmessage += std::to_string(prtId) + " which has PDG: " + std::to_string(prtPdg);
  Log(logmessage, v_message, verbosity);

}

void BackTracker::DirectParentsFromClockTickWindows()
{
  if (!fHitToDirectParents) return;
  fHitToDirectParents->clear();
  if (!fPMTToDirectParentMap || !fRecoADCHits) return;

  const double prewindow_ns = static_cast<double>(fPMTSimPrewindowTicks) * NS_PER_ADC_SAMPLE;
  const double readout_ns = static_cast<double>(fPMTSimReadoutWindowTicks) * NS_PER_ADC_SAMPLE;

  uint32_t evtNum = 0;
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber", evtNum);
  std::cout << "BackTracker::DirectParentsFromClockTickWindows [event=" << evtNum << "] input:"
            << " RecoADCHits PMTs=" << fRecoADCHits->size()
            << ", PMTToDirectParentMap PMTs=" << fPMTToDirectParentMap->size()
            << ", preTicks=" << fPMTSimPrewindowTicks
            << ", readoutTicks=" << fPMTSimReadoutWindowTicks
            << std::endl;

  size_t total_pulses_checked = 0;
  size_t total_tick_matches = 0;
  size_t total_parent_ids_added = 0;
  std::set<int> unique_parent_ids;
  int debug_print_budget = 20;

  for (auto const& recoIt : *fRecoADCHits) {
    unsigned long pmtID = recoIt.first;

    auto parentMapIt = fPMTToDirectParentMap->find(pmtID);
    if (parentMapIt == fPMTToDirectParentMap->end()) continue;

    std::map<uint16_t, std::vector<int>> const& hits_to_directparents_map = parentMapIt->second;

    std::cout << "BackTracker::DirectParentsFromClockTickWindows PMT " << pmtID
              << ": reco minibufs=" << recoIt.second.size()
              << ", truth tick bins=" << hits_to_directparents_map.size()
              << std::endl;

    for (std::vector<ADCPulse> const& minibufPulses : recoIt.second) {
      for (ADCPulse const& pulse : minibufPulses) {
        ++total_pulses_checked;
        double pulseStart = pulse.start_time();
        double hitTime = pulse.peak_time();

        double tmin = pulseStart - prewindow_ns;
        double tmax = pulseStart + readout_ns;

        if (debug_print_budget > 0) {
          std::cout << "BackTracker::DirectParentsFromClockTickWindows pulse: PMT=" << pmtID
                    << ", start=" << pulseStart
                    << ", peak=" << hitTime
                    << ", tmin=" << tmin
                    << ", tmax=" << tmax
                    << std::endl;
          --debug_print_budget;
        }

        for (auto const& apair : hits_to_directparents_map) {
          double mchitTime = static_cast<double>(apair.first) * NS_PER_ADC_SAMPLE;

          if (mchitTime > tmin && mchitTime < tmax) {
            ++total_tick_matches;
            (*fHitToDirectParents)[pmtID][hitTime].insert(
              (*fHitToDirectParents)[pmtID][hitTime].end(),
              apair.second.begin(), apair.second.end());
            total_parent_ids_added += apair.second.size();
            unique_parent_ids.insert(apair.second.begin(), apair.second.end());

            if (debug_print_budget > 0) {
              std::cout << "BackTracker::DirectParentsFromClockTickWindows match: PMT=" << pmtID
                        << ", mchitTime=" << mchitTime
                        << ", tick=" << apair.first
                        << ", nParentsAdded=" << apair.second.size()
                        << std::endl;
              --debug_print_budget;
            }
          }
        }
      }
    }
  }

  std::cout << "BackTracker::DirectParentsFromClockTickWindows output: matched PMTs="
            << fHitToDirectParents->size()
            << ", pulsesChecked=" << total_pulses_checked
            << ", tickMatches=" << total_tick_matches
            << ", parentIDsAdded=" << total_parent_ids_added
            << ", uniqueParentIDs=" << unique_parent_ids.size()
            << std::endl;

  if (!fHitToDirectParents->empty()) {
    auto it = fHitToDirectParents->begin();
    std::cout << "BackTracker::DirectParentsFromClockTickWindows example: PMT " << it->first
              << " has " << it->second.size() << " hit times with direct parent matches."
              << std::endl;
  }
}


//------------------------------------------------------------------------------
bool BackTracker::LoadFromStores()
{
  // Grab the stuff we need from the stores
  bool goodMCClusters = m_data->CStore.Get("ClusterMapMC", fClusterMapMC);
  if (!goodMCClusters) {
    std::cerr<<"BackTracker: no ClusterMapMC in the CStore!"<<endl;
    return false;
  }

  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    std::cerr<<"BackTracker: no ANNIEEvent store!"<<endl;
    return false;
  }
    
  bool goodMCHits = m_data->Stores.at("ANNIEEvent")->Get("MCHits", fMCHitsMap);
  if (!goodMCHits) {
    std::cerr<<"BackTracker: no MCHits in the ANNIEEvent!"<<endl;
    return false;
  }
  
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    std::cerr<<"BackTracker: no MCParticles in the ANNIEEvent!"<<endl;
    return false;
  }

  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    std::cerr<<"BackTracker: no TrackId_to_MCParticleIndex in the ANNIEEvent!"<<endl;
    return false;
  }

  if (fUsePulseWindowMatching) {
    bool gotDirectParentMap = m_data->Stores.at("ANNIEEvent")->Get("PMTToDirectParentMap", fPMTToDirectParentMap);
    bool gotRecoADCHits = m_data->Stores.at("ANNIEEvent")->Get("RecoADCHits", fRecoADCHits);
    if (!gotDirectParentMap || !gotRecoADCHits) {
      logmessage = "BackTracker: PMTToDirectParentMap or RecoADCHits missing, disabling pulse-window matching for this event.";
      Log(logmessage, v_warning, verbosity);
      fPMTToDirectParentMap = nullptr;
      fRecoADCHits = nullptr;
    }

    uint16_t prewindowTicks = fPMTSimPrewindowTicks;
    uint16_t readoutTicks = fPMTSimReadoutWindowTicks;
    if (m_data->Stores.at("ANNIEEvent")->Get("PMTSimPrewindowTicks", prewindowTicks)) {
      fPMTSimPrewindowTicks = prewindowTicks;
    }
    if (m_data->Stores.at("ANNIEEvent")->Get("PMTSimReadoutWindowTicks", readoutTicks)) {
      fPMTSimReadoutWindowTicks = readoutTicks;
    }
  }

  return true;
}
