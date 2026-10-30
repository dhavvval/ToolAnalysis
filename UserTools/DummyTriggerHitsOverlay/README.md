# DummyTriggerHitsOverlay

DummyTriggerHitsOverlay adds the tank hits of one dummy-trigger (random-trigger) data window to every simulated event, so that MC events carry the same accidental tank activity as data (dark noise, radioactivity, afterpulses, cosmic remnants) when they are clustered.

The hits are added to the waveform MC hits, after `PhaseIIADCHitFinder` (`BuildWaveformMCHits 1`) and before `ClusterFinder`. Clustering, the cluster time merge and the truth matching then see the simulated and the data hits together, as for a data trigger.

## Tool order

```
LoadWCSim -> PMTWaveformSim -> PhaseIIADCHitFinder (MCWaveforms 1, BuildWaveformMCHits 1)
          -> DummyTriggerHitsOverlay -> ClusterFinder (HitStore WaveformMCHits) -> BackTracker -> ANNIEEventTreeMaker
```

## Input

`OverlayFile` is a BeamClusterAnalysis ntuple (tree `Event`) of a dummy-trigger run, e.g. `BeamCluster_6256.root` (LAPPD off) or `BeamCluster_6254.root` (LAPPD on). The tool reads `eventTimeTank`, `hitT`, `hitPE` and `hitChankey`, which hold every tank hit of the trigger window.

As in `ambe data randoms`, only triggers with tank data are used (`eventTimeTank` != 0), one per `eventTimeTank`.

## What is added

For each MC event the tool picks one dummy trigger. The pick is a hash of `Seed`, the WCSim file name (without path) and the WCSim entry number, so it is reproducible and does not depend on job splitting or event order. Triggers are drawn with replacement. The file name is the one LoadWCSim reads: the AmBe grid jobs link each file as `wcsim_0.<runNumber>.<file index>.root`, and runNumber carries the source position, so different positions draw different triggers.

Each data hit becomes one `MCHit` in `WaveformMCHits`:
- channel: data channel key -> WCSim tube ID (`channelkey_to_pmtid_data`) -> MC channel key (`pmt_tubeid_to_channelkey`); hits on channels without an MC channel or on channels that are OFF are dropped;
- time: `hitT + TimeOffset` [ns];
- charge: `hitPE` [p.e.], the unit of the waveform MC hits (pulse charge / SPE gain, same gain table);
- `Parents` empty, `DirectParents = {-3}` (overlay; -1 is dark noise, -2 a track WCSim did not save), `IsDarknoise` false.

`WaveformMCHitNTrue` and `WaveformMCHitTrueCharge` stay aligned with the hit store: overlay hits get 0 true MCHits and 0 true p.e., so `hitTruePE` is 0 for them. Overlay hits are appended at the end of each channel's list, so `MCHitToWaveformHit` indices stay valid.

The ANNIEEvent also gets `DummyTriggerOverlayEntry` (tree entry of the trigger used), `DummyTriggerOverlayTankTime` and `DummyTriggerOverlayNHits`.

Hits are added as hits, not as waveforms: two pulses on the same PMT within a few ns stay two hits. With about 44 accidental hits per 70 us window over ~120 PMTs (6256) this pile-up is negligible.

## Time frame

MC hit times start at the primary vertex (t = 0 is the neutron emission). A data window starts before its trigger: in AmBe source runs the prompt tank light sits 1,094–1,104 ns into the window (median hit time, 12 runs March–July 2026, all ports), while the MC prompt light sits at 23–25 ns, so the two clocks differ by about 1,074 ns.

`TimeOffset` must match the clock the MC hits are in:
- MC in its own clock (PMTWaveformSim `TimeShift 0`): `TimeOffset -1074` puts the dummy-window hits where they sit relative to the emission in a real AmBe trigger.
- MC moved to the data clock (PMTWaveformSim `TimeShift 1074`): `TimeOffset 0`. ClusterFinder's end-of-window cut and the time-merge threshold then act at the same window times as in data. This is the AmBe production setting (grid_wcsim `--chain wf_dtho`).

The value depends on the trigger and DAQ settings of the data and on the MC time convention, so measure it again for other triggers (beam, cosmic, LED/laser) or other MC set-ups.

## Configuration

```
verbosity 1
OverlayFile ./BeamCluster_6256.root   # dummy-trigger BeamCluster ntuple
TreeName Event                        # default Event
HitStore WaveformMCHits               # only WaveformMCHits is supported
Seed 0                                # pick = hash(Seed, WCSim file name, WCSim entry)
TimeOffset -1074                      # ns, added to the data hit times; -1074 with PMTWaveformSim TimeShift 0, 0 with TimeShift 1074 (AmBe)
```
