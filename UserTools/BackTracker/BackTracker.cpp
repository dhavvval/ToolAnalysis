#include "BackTracker.h"

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


  // Set up the pointers we're going to save. No need to 
  // delete them at Finalize, the store will handle it
  fClusterToBestParticleID  = new std::map<double, int>;
  fClusterToBestParticlePDG = new std::map<double, int>;
  fClusterEfficiency        = new std::map<double, double>;
  fClusterPurity            = new std::map<double, double>;
  fClusterTotalCharge       = new std::map<double, double>;
  fClusterHitToDirectParentTrackIDs = new std::map<double, std::vector<std::vector<int>>>;
  
  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Execute()
{
  if (!LoadFromStores())
    return false;

  fClusterToBestParticleID ->clear();
  fClusterToBestParticlePDG->clear();
  fClusterEfficiency       ->clear();
  fClusterPurity           ->clear();
  fClusterTotalCharge      ->clear();

  fParticleToTankTotalCharge.clear();
  fClusterHitToDirectParentTrackIDs->clear();
  SumParticleTankCharge();

  //mapping MC-Index to TrackIDs
  std::map<int,int> MCIndexToTrackID;
  for (const auto& mchit : *fMCParticleIndexMap) {
    MCIndexToTrackID[mchit.second] = mchit.first; 
  }

  std::cout << "[BT DEBUG] MCIndexToTrackID size = " << MCIndexToTrackID.size() << std::endl;

  // Loop over the clusters and do the things
  for (std::pair<double, std::vector<MCHit>>&& apair : *fClusterMapMC) {
    int prtId = -5;
    int prtPdg = -5;
    double eff = -5;
    double pur = -5;
    double totalCharge = 0;

    // Grabbing DirectparentIDs for each MCHits in the cluster
    std::cout << "[BT DEBUG] cluster_time=" << apair.first
          << " mchits_in_cluster=" << apair.second.size() << std::endl;

    std::vector<std::vector<int>> clusterMCHits_DirectIDs;
    
    clusterMCHits_DirectIDs.reserve(apair.second.size());
    for (auto& mchit : apair.second) {
      std::vector<int> directParentIDs;
      const std::vector<int>* directIdxs = mchit.GetDirectParents();

      for (int idx : *directIdxs) {
        auto it = MCIndexToTrackID.find(idx);
        if (it != MCIndexToTrackID.end()) directParentIDs.push_back(it->second);
      }

      std::cout << "[BT DEBUG] hit directIdx_count=" << directIdxs->size()
          << " mapped_directTrackID_count=" << directParentIDs.size() << " IDs: ";
          for (int id : directParentIDs) std::cout << id << " ";
          std::cout << std::endl;

      clusterMCHits_DirectIDs.push_back(std::move(directParentIDs));
    }

    MatchMCParticle(apair.second, prtId, prtPdg, eff, pur, totalCharge);

    fClusterToBestParticleID ->emplace(apair.first, prtId);
    fClusterToBestParticlePDG->emplace(apair.first, prtPdg);
    fClusterEfficiency       ->emplace(apair.first, eff);
    fClusterPurity           ->emplace(apair.first, pur);
    fClusterTotalCharge      ->emplace(apair.first, totalCharge);
    fClusterHitToDirectParentTrackIDs->emplace(apair.first, std::move(clusterMCHits_DirectIDs));
    std::cout << "[BT DEBUG] stored vectors for cluster_time=" << apair.first
          << " count=" << fClusterHitToDirectParentTrackIDs->at(apair.first).size()
          << std::endl;

  }

  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticleID",  fClusterToBestParticleID );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  m_data->Stores.at("ANNIEEvent")->Set("ClusterEfficiency",        fClusterEfficiency       );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterPurity",            fClusterPurity           );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterTotalCharge",       fClusterTotalCharge      );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterHitToDirectParentTrackIDs", fClusterHitToDirectParentTrackIDs);

std::map<double, std::vector<std::vector<int>>>* check = nullptr;
bool ok = m_data->Stores.at("ANNIEEvent")->Get("ClusterHitToDirectParentTrackIDs", check);
std::cout << "[BT DEBUG] ANNIEEvent Get ClusterHitToDirectParentTrackIDs ok=" << ok
          << " clusters=" << (ok ? check->size() : 0) << std::endl;

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

  return true;
}
