#include "DummyTriggerHitsOverlay.h"

#include <set>
#include <sstream>

#include "TFile.h"
#include "TTree.h"

namespace {

uint64_t SplitMix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t Fnv1a(const std::string& s) {
  uint64_t h = 1469598103934665603ULL;
  for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
  return h;
}

}

DummyTriggerHitsOverlay::DummyTriggerHitsOverlay():Tool(){}


bool DummyTriggerHitsOverlay::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  m_data= &data;

  m_variables.Get("verbosity",verbosity);
  if (!m_variables.Get("OverlayFile",fOverlayFile)) {
    Log("DummyTriggerHitsOverlay: OverlayFile not set. Aborting!", v_error, verbosity);
    return false;
  }
  m_variables.Get("TreeName",fTreeName);
  m_variables.Get("HitStore",fHitStoreName);
  m_variables.Get("Seed",fSeed);
  m_variables.Get("TimeOffset",fTimeOffset);

  if (fHitStoreName != "WaveformMCHits") {
    Log("DummyTriggerHitsOverlay: only HitStore WaveformMCHits is supported. Aborting!", v_error, verbosity);
    return false;
  }
  if (!m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",geom)) {
    Log("DummyTriggerHitsOverlay: no AnnieGeometry in the ANNIEEvent header (run LoadWCSim first). Aborting!", v_error, verbosity);
    return false;
  }
  if (!LoadTriggers() || !BuildChannelMap()) return false;
  fTimesUsed.assign(fTrigEntry.size(), 0);

  return true;
}


bool DummyTriggerHitsOverlay::LoadTriggers(){

  TFile* f = TFile::Open(fOverlayFile.c_str(), "READ");
  if (f == nullptr || f->IsZombie()) {
    Log("DummyTriggerHitsOverlay: cannot open " + fOverlayFile + ". Aborting!", v_error, verbosity);
    return false;
  }
  TTree* t = (TTree*) f->Get(fTreeName.c_str());
  const char* names[] = {"eventTimeTank", "hitT", "hitPE", "hitChankey"};
  bool ok = (t != nullptr);
  for (const char* b : names) ok = ok && t->GetBranch(b) != nullptr;
  if (!ok) {
    Log("DummyTriggerHitsOverlay: " + fOverlayFile + " has no tree " + fTreeName + " with eventTimeTank, hitT, hitPE, hitChankey. Aborting!", v_error, verbosity);
    f->Close();
    delete f;
    return false;
  }

  ULong64_t ett = 0;
  std::vector<double>* ht = nullptr;
  std::vector<double>* hpe = nullptr;
  std::vector<int>* hck = nullptr;
  t->SetBranchStatus("*", 0);
  for (const char* b : names) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("eventTimeTank", &ett);
  t->SetBranchAddress("hitT", &ht);
  t->SetBranchAddress("hitPE", &hpe);
  t->SetBranchAddress("hitChankey", &hck);

  // as `ambe data randoms`: triggers with tank data only, one per eventTimeTank
  std::set<uint64_t> seen;
  long nentries = t->GetEntries();
  fTrigFirstHit.push_back(0);
  for (long i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    if (ett == 0 || !seen.insert(ett).second) continue;
    for (size_t k = 0; k < ht->size(); ++k) {
      fHitChankey.push_back(hck->at(k));
      fHitT.push_back(ht->at(k));
      fHitPE.push_back(hpe->at(k));
    }
    fTrigEntry.push_back(i);
    fTrigTankTime.push_back(ett);
    fTrigFirstHit.push_back(fHitT.size());
  }
  f->Close();
  delete f;

  if (fTrigEntry.empty()) {
    Log("DummyTriggerHitsOverlay: no trigger with tank data in " + fOverlayFile + ". Aborting!", v_error, verbosity);
    return false;
  }
  std::stringstream ss;
  ss << "DummyTriggerHitsOverlay: " << fTrigEntry.size() << " dummy triggers (" << nentries << " entries) and "
     << fHitT.size() << " hits from " << fOverlayFile;
  Log(ss.str(), v_message, verbosity);
  return true;
}


bool DummyTriggerHitsOverlay::BuildChannelMap(){

  std::map<unsigned long, int> data_chankey_to_tubeid;
  std::map<int, unsigned long> tubeid_to_mc_chankey;
  if (!m_data->CStore.Get("channelkey_to_pmtid_data", data_chankey_to_tubeid)
      || !m_data->CStore.Get("pmt_tubeid_to_channelkey", tubeid_to_mc_chankey)) {
    Log("DummyTriggerHitsOverlay: channelkey_to_pmtid_data / pmt_tubeid_to_channelkey not in the CStore (run LoadWCSim first). Aborting!", v_error, verbosity);
    return false;
  }
  for (auto const& apair : data_chankey_to_tubeid) {
    auto it = tubeid_to_mc_chankey.find(apair.second);
    if (it != tubeid_to_mc_chankey.end()) fDataToMCChankey[apair.first] = it->second;
  }
  return true;
}


size_t DummyTriggerHitsOverlay::PickTrigger(const std::string& mcfile, uint64_t mcevent) const {
  std::string base = mcfile.substr(mcfile.find_last_of('/') + 1);
  uint64_t key = SplitMix64(SplitMix64(Fnv1a(base) ^ (uint64_t)(int64_t)fSeed) ^ mcevent);
  return key % fTrigEntry.size();
}


bool DummyTriggerHitsOverlay::Execute(){

  BoostStore* annie_event = m_data->Stores["ANNIEEvent"];
  std::string mcfile;
  uint64_t mcevent = 0;
  if (!annie_event->Get("MCFile", mcfile) || !annie_event->Get("MCEventNum", mcevent)) {
    Log("DummyTriggerHitsOverlay: no MCFile / MCEventNum in the ANNIEEvent. Aborting!", v_error, verbosity);
    return false;
  }
  std::map<unsigned long, std::vector<MCHit>>* wf_mchits = nullptr;
  if (!annie_event->Get(fHitStoreName, wf_mchits)) {
    Log("DummyTriggerHitsOverlay: no " + fHitStoreName + " (run PhaseIIADCHitFinder with BuildWaveformMCHits 1 first). Aborting!", v_error, verbosity);
    return false;
  }
  // kept aligned with the hit store: overlay hits have no true MCHits
  std::map<unsigned long, std::vector<int>> wf_true_n;
  std::map<unsigned long, std::vector<double>> wf_true_charge;
  annie_event->Get("WaveformMCHitNTrue", wf_true_n);
  annie_event->Get("WaveformMCHitTrueCharge", wf_true_charge);

  size_t itrig = PickTrigger(mcfile, mcevent);
  int n_added = 0;
  for (size_t k = fTrigFirstHit[itrig]; k < fTrigFirstHit[itrig + 1]; ++k) {
    auto it = fDataToMCChankey.find(fHitChankey[k]);
    if (it == fDataToMCChankey.end()) { ++fNUnmapped; continue; }
    unsigned long chankey = it->second;
    Channel* thischannel = geom->GetChannel(chankey);
    if (thischannel == nullptr || thischannel->GetStatus() == channelstatus::OFF) { ++fNChannelOff; continue; }

    std::vector<MCHit>& hits = (*wf_mchits)[chankey];
    std::vector<int>& truen = wf_true_n[chankey];
    std::vector<double>& trueq = wf_true_charge[chankey];
    truen.resize(hits.size(), 0);
    trueq.resize(hits.size(), 0.);
    hits.push_back(MCHit(chankey, fHitT[k] + fTimeOffset, fHitPE[k], std::vector<int>{}, std::vector<int>{kOverlayDirectParent}));
    truen.push_back(0);
    trueq.push_back(0.);
    ++n_added;
  }
  annie_event->Set("WaveformMCHitNTrue", wf_true_n);
  annie_event->Set("WaveformMCHitTrueCharge", wf_true_charge);
  annie_event->Set("DummyTriggerOverlayEntry", (int) fTrigEntry[itrig]);
  annie_event->Set("DummyTriggerOverlayTankTime", fTrigTankTime[itrig]);
  annie_event->Set("DummyTriggerOverlayNHits", n_added);

  ++fNEvents;
  fNAdded += n_added;
  ++fTimesUsed[itrig];
  if (verbosity >= v_debug) {
    std::stringstream ss;
    ss << "DummyTriggerHitsOverlay: " << mcfile << " event " << mcevent << " <- dummy trigger entry "
       << fTrigEntry[itrig] << " (" << n_added << " hits)";
    Log(ss.str(), v_debug, verbosity);
  }
  return true;
}


bool DummyTriggerHitsOverlay::Finalise(){

  int n_distinct = 0;
  int max_used = 0;
  for (int n : fTimesUsed) {
    if (n > 0) ++n_distinct;
    if (n > max_used) max_used = n;
  }
  std::stringstream ss;
  ss << "DummyTriggerHitsOverlay summary: " << fNEvents << " MC events, " << fNAdded << " overlay hits ("
     << (fNEvents ? double(fNAdded) / fNEvents : 0.) << " per event), " << n_distinct << " of "
     << fTrigEntry.size() << " dummy triggers used (at most " << max_used << " times), hits dropped: "
     << fNUnmapped << " no MC channel, " << fNChannelOff << " channel off";
  Log(ss.str(), v_warning, verbosity);
  return true;
}
