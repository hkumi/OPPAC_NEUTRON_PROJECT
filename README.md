# Geant4 Simulation of a Gas Scintillation Neutron Imaging Detector with 1 mm SiPM Sensors

> **Phase 1** — Proof-of-concept simulation of a novel neutron imaging detector  
> using Ar:CF4 gas scintillation and 1 mm pitch Silicon Photomultiplier (SiPM) arrays.

---

## Overview

This repository contains a Geant4-based Monte Carlo simulation of a fast neutron imaging detector designed for medical device inspection. The detector exploits gas scintillation in an Ar:CF4 (90:10) mixture at 30 mbar, with optical photon readout via four radial arrays of 1 mm pitch SiPMs.

To our knowledge, **this is the first simulation study to investigate 1 mm pitch SiPM sensors for gas scintillation neutron imaging**, establishing a spatial resolution baseline of **0.52 mm** at 2.5 MeV neutron energy.

A valve phantom (stainless steel body with water cavity and aluminium stem) is used as the imaging target, representative of internal medical device components that are inaccessible to conventional X-ray imaging.

---

## Key Features

- 2.5 MeV neutron beam (D-D fusion source energy)
- Ar:CF4 90:10 scintillating gas at 30 mbar
- HDPE neutron-to-proton recoil converter (1.5 mm)
- Four radial SiPM arrays — 25 sensors per side, **1 mm pitch** (novel)
- Teflon collimator with optical reflectivity 95%
- Mylar mirror electrodes for optical confinement
- Gaussian beam profile (σ = 5 mm)
- Full optical photon tracking with realistic SiPM quantum efficiency (25–50%)
- Angle-dependent spatial resolution analysis
- Gaussian MLE centroid reconstruction algorithm
- Valve phantom: steel body, water cavity, aluminium stem

---

## Physics

### Detection Principle

```
2.5 MeV neutron
      ↓
HDPE converter (n,p elastic scattering)
      ↓
Proton recoil enters Ar:CF4 gas
      ↓
Gas scintillation (UV/visible photons)
      ↓
Optical photons propagate to radial SiPM arrays
      ↓
Centroid reconstruction → 2D interaction position
      ↓
Neutron image of phantom
```

### Gas Properties

| Parameter | Value |
|-----------|-------|
| Mixture | Ar:CF4 90:10 |
| Pressure | 30 mbar |
| Scintillation yield | 2.5 × 10⁶ photons/MeV |
| Time constant | 15 ns |

### Spatial Resolution

| Metric | Value |
|--------|-------|
| Spatial resolution (FWHM) | **0.52 mm** |
| SiPM pitch | 1.0 mm |
| Reconstruction method | Weighted Gaussian MLE |
| Neutron energy | 2.5 MeV |

---

## Detector Geometry

```
z = -20 cm   z = -15 cm        z = -0.5 cm    z = 0     z = +0.5 cm
  [source] →  [valve phantom] →  [gas box entrance]       [gas box exit]
                                  [MylK cathode]            [MylA anode]
                                  [      Ar:CF4 gas + Teflon collimator      ]
                                  [HDPE converter at z=+1.3mm inside gas]
                                  [SiPMs: radial at z=0, 4×25 = 100 sensors]
```

### Volumes

| Volume | Material | Dimensions |
|--------|----------|------------|
| World | Vacuum | 60×60×100 cm |
| Detector gas box | Ar:CF4 30 mbar | 10×10×1 cm |
| HDPE converter | Polyethylene | 27.5×27.5×1.5 mm |
| Collimator cells | Teflon | 3.0×1.1×1.1 mm (×100) |
| SiPM sensors | SiO₂ approx. | 1.0×1.0×1.0 mm (×100) |
| Mylar electrodes | Aluminium | 27.5×27.5×0.012 mm |
| Valve body | Stainless steel | R=8 mm, L=4 mm |
| Valve cavity | Water | R=4.99 mm, L=3.98 mm |
| Valve stem | Aluminium | R=1.95 mm, L=3.5 mm |

---

## Repository Structure

```
├── src/
│   ├── DetectorConstruction.cc   # Geometry and materials
│   ├── PrimaryGeneratorAction.cc # 2.5 MeV neutron beam
│   ├── SteppingAction.cc         # True position capture + energy tracking
│   ├── EventAction.cc            # Per-event data management
│   ├── RunAction.cc              # Histogram/ntuple booking
│   └── detector.cc               # SiPM sensitive detector + reconstruction
├── include/
│   ├── DetectorConstruction.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── SteppingAction.hh
│   ├── EventAction.hh
│   ├── RunAction.hh
│   ├── detector.hh
│   └── SensorHit.hh
├── CMakeLists.txt
├── run.mac                        # Batch run macro
├── vis.mac                        # Visualisation macro
└── README.md
```

---

## Requirements

| Dependency | Version |
|------------|---------|
| Geant4 | ≥ 11.0 (with optical physics) |
| CMake | ≥ 3.16 |
| C++ | ≥ C++17 |
| ROOT | ≥ 6.24 (for analysis output) |

Geant4 must be built with the following options enabled:
```
GEANT4_USE_OPENGL_X11=ON   (for visualisation)
GEANT4_BUILD_MULTITHREADED=ON
```

---

## Build Instructions

```bash
# Clone the repository
git clone (https://github.com/hkumi/OPPAC_NEUTRON_PROJECT.git)
cd OPPAC_NEUTRON_PROJECT

# Create build directory
mkdir build && cd build

# Configure
cmake .. -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4

# Build
make -j$(nproc)
```

---

## Running the Simulation

### Batch mode
```bash
./example4a -m run.mac
```

### Interactive / Visualisation mode
```bash
./example4a
```

### Recommended run.mac
```
/run/numberOfThreads 4
/run/initialize
/run/printProgress 1000
/run/beamOn 100000
```

---

## Output

The simulation produces a ROOT file containing:

| Object | ID | Description |
|--------|----|-------------|
| H1 | 0 | Optical photon energy spectrum |
| H1 | 1–4 | SiPM hit distributions (top/bottom/left/right) |
| H1 | 6–7 | Reconstructed X and Y position |
| H1 | 8–9 | Position residuals ΔX, ΔY |
| H1 | 10–12 | Proton energy in converter/cathode/gas |
| H1 | 13 | Spatial resolution distribution |
| H2 | 0 | Reconstructed 2D image |
| H2 | 2–4 | Theta vs ΔX, ΔY, resolution |
| H2 | 5–6 | Phi vs ΔX, ΔY |
| Ntuple | 5 | Full per-event reconstruction data |

### Key ntuple columns (ID 5)

| Column | Description |
|--------|-------------|
| x_rec, y_rec | Reconstructed position (mm) |
| x_true, y_true | True proton entry position (mm) |
| dx, dy | Position residuals (mm) |
| theta, phi | Proton direction angles (degrees) |
| resolution | 2D spatial resolution (mm) |

---

## Physics List

Recommended physics list for 2.5 MeV neutron imaging:

```
/physics/addPhysics FTFP_BERT_HP
```

The `_HP` suffix enables the **high-precision neutron package** (< 20 MeV), essential for accurate neutron transport and elastic scattering cross-sections.

Add optical physics for scintillation and photon tracking:
```
/physics/addPhysics G4OpticalPhysics
```

---

## Reconstruction Algorithm

Position reconstruction uses a **Weighted Gaussian Maximum Likelihood Estimator (MLE)**:

1. For each SiPM array, compute centroid from hit sensor indices
2. Weight each array's estimate by N/σ (photon count / spread)
3. Combine opposing arrays (top+bottom for X, left+right for Y)

```
x_rec = (mean_top × N_top/σ_top + mean_bottom × N_bottom/σ_bottom)
        / (N_top/σ_top + N_bottom/σ_bottom)
```

Dark noise is simulated as a Poisson process (λ = 6 counts/event/array).

---

## Results Summary (Phase 1)

| Parameter | Value |
|-----------|-------|
| Neutron energy | 2.5 MeV |
| Spatial resolution | **0.52 mm** |
| SiPM pitch | 1 mm (novel — first study at this pitch) |
| Gas pressure | 30 mbar |
| Converter thickness | 1.5 mm HDPE |
| Beam profile | Gaussian, σ = 5 mm |
| Phantom | Medical valve (steel/water/aluminium) |

---

## Novelty and Publication Plan

### Phase 1 (This Repository)
**Novel contribution:** First Monte Carlo study of 1 mm pitch SiPM sensors for Ar:CF4 gas scintillation neutron imaging.

- Establishes 0.52 mm spatial resolution baseline
- Demonstrates valve phantom imaging capability
- Validates reconstruction algorithm for sub-mm SiPM pitch

### Phase 2 (will be done in future)
**Planned contribution:** Systematic optimisation study including:
- SiPM pitch scan (0.5 mm – 3 mm)
- Gas pressure optimisation (30 mbar -300mbar)
- Converter thickness study (1 mm – 10 mm HDPE)
- Angle-dependent resolution characterisation
- Comparison with experimental measurements

---


---

## Acknowledgements

This work was carried out as part of a PhD research project in neutron detector development.  
Simulation framework based on the Geant4 B4 example.

---

## License

The Geant4 simulation framework is used under the [Geant4 Software License](http://cern.ch/geant4/license).

---

## Contact

For questions or collaboration enquiries, please open an issue or contact:  
**Harriet Kumi** — h.kumi@udc.es
