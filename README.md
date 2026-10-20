# JUNO // Event Display

An interactive 3D event display for the JUNO detector that runs entirely in a web browser.

- Built-in generator for atmospheric νμ / νe CC, atmospheric NC, reactor IBD and cosmic muons, with tunable energy, direction, flavour, vertex and number of events.
- Custom mode: aim the neutrino (zenith, azimuth), set its energy, and pick the interaction (QE, MEC/2p2h, RES, DIS, coherent π, or automatic by energy).
- Opens GENIE `.gst.root` files directly in the browser (no server, no ROOT install).
- Optional Honda (HAKKM) flux file, the same one used by `gevgen_atmo`.
- Simplified light model: scintillation + Cherenkov light, 17,596 20" PMTs, 25,587 3" PMTs, 3,352 water-pool PMTs, multiplicity trigger, approximate vertex/energy/track reconstruction.
- Export of truth particles, PMT hits (JSON/CSV) and per-run summaries for analysis.

## Run it

Open the GitHub Pages link of this repository. No installation is needed.

To run it on your own computer instead, the page must be served (opening `index.html` by double-click does not work because of browser security rules for modules):

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

## Use your GENIE files

1. Click **GENIE file** in the Generator panel.
2. Click **Open .gst.root file** (or drag the file onto the 3D view).
3. Choose the first entry and how many to read (up to 1000).

The file never leaves your computer. Vertices in gst are in metres (SI) and are used as JUNO coordinates (origin at the CD centre, z up). If the vertices look like centimetres, the page converts them automatically.

## Keyboard

| Key | Action |
|---|---|
| Space | Play / pause |
| N | Next event |
| R | Reset camera |
| 1–6 | Camera presets |
| C | Cycle geometry mode |
| L | Labels on/off |
| Esc | Clear selection |

## Exported data

- **Event JSON** (`juno-event-display/1`): event header, all truth particles (PDG, 4-momentum, start, end, fate, energy deposited and visible in LS), and every readout window with all PMT hits (`pmt`, `type`, `npe`, `t_ns`, `dark`) and the approximate reconstruction.
- **Hits CSV**: one row per fired PMT with θ, φ, position, p.e. and first-hit time.
- **Run CSV / JSON**: one row per event in the run with neutrino, interaction, vertex, deposited and visible energy, and an expected prompt p.e. estimate.

Units: metres, nanoseconds, GeV (fields ending in `_MeV` are MeV).

## Geometry sources

Detector values come from:

- JUNO Collaboration, *First measurement of reactor neutrino oscillations at JUNO*, arXiv:2511.14593.
- JUNO Collaboration, *Initial performance results of the JUNO detector*, arXiv:2511.14590.

Used directly: AV inner radius 17.7 m (263 panels, 120 mm, 590 bars), SS structure inner diameter 40.1 m, 17,596 20" PMTs (4,939 dynode, 12,657 MCP), 25,587 3" PMTs, water pool Ø 43.5 m × 44 m, 2,404 + 348 20" and 600 8" water-pool PMTs (ring heights from the paper), Top Tracker 3 layers × 7 × 3 walls of 6.8 m, 1785 p.e./MeV at the centre, LS attenuation 20.6 m, TTS 1.93 / 7.53 ns, dark rates 20.6 / 22.7 kHz, 1008 ns gate, trigger 350 PMTs in 304 ns, fiducial cut R < 16.5 m, |z| < 15.5 m, effective refractive index 1.60.

All geometry constants are in `js/detector.js`, each marked `[P]` (from the papers) or `[D]` (display approximation).

## Scientific approximations and limitations

This is an educational display, not JUNOSW, Geant4 or real JUNO data.

- **PMT positions** are a deterministic quasi-uniform layout with the correct counts and coverage, not the engineering map. Vertical placement of the pool, TT heights and chimney are approximate.
- **Built-in generator**: free-nucleon kinematics (no Fermi motion, no FSI) with simple QE / MEC (nucleon pair) / Δ resonance / DIS / coherent-π models; the automatic process mix uses approximate CC cross-section shapes, not a GENIE tune; cross sections are not modelled, so event rates and process fractions are illustrative. Extra low-energy "nuclear" neutrons are added for illustration and are flagged in the particle list. GENIE files give the real generator output.
- **Particle transport**: straight tracks with continuous energy loss (no magnetic field in JUNO), parameterised EM showers (X₀ = 0.52 m), Compton chains for MeV γ, simple hadronic interactions, stopped μ⁻ decay (τ = 2.026 µs) or are captured on ¹²C (7.8%, emitting neutrons), stopped μ⁺ decay (τ = 2.197 µs), a neutron random walk ending in capture on H (2.223 MeV) or ¹²C (4.946 MeV) after τ ≈ 200 µs. Not a substitute for Geant4.
- **Light**: isotropic scintillation with Birks-like quenching for protons and hadronic deposits, a 6% directional Cherenkov component on a 48° cone, straight-line propagation with attenuation, no refraction or reflection at the acrylic, no re-emission. LS energy non-linearity is not modelled.
- **Reconstruction**: charge-weighted vertex with a bias correction, energy from total p.e., and a two-cluster entry/exit fit for muon-like events. These are simple estimates to compare with truth, not JUNO's OMILREC / VTREP / JVertex.
- **Cosmic muons**: approximate underground spectrum (mean ~200 GeV) and zenith distribution; the mountain profile is not included.
- **Timing**: animation is slowed educational playback. Delayed windows (Michel decay, neutron captures) are shown with compressed Δt.
- Downloads are blocked inside some embedded previews; they work on the GitHub Pages site.

## Repository layout

The live site is the single `index.html` at the top of the repository; it contains all the code. The `source/` folder holds the same code as separate, readable files for editing. After editing them, rebuild `index.html` (ask Claude, or run `python3 build_single.py` with esbuild installed).

## Source files

```
source/index.html  page layout and styles
js/app.js         interface, animation, rendering loop
js/scene.js       Three.js detector model and event visuals
js/physics.js     generators, propagation, light/PMT response, reconstruction
js/detector.js    geometry and detector constants
js/plots.js       PMT map, time histogram, projections
js/gst.js         GENIE gst reader (JSROOT)
js/worker.js      background thread for the PMT response
```

Libraries load from cdn.jsdelivr.net: three.js 0.160.0 and JSROOT 7.10.3.
