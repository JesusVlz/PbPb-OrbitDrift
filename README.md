# PbPb Orbit Drift

Notebook-based analysis of the residual orbit drift (ROD) during Pbâ€“Pb van der Meer scans, using CMS DOROS and arc BPM measurements. The current configuration and example inputs are for LHC fill **9212** (2023, Pbâ€“Pb at 5.36 TeV per nucleon pair).

The workflow reads nominal scan positions and BPM time series, estimates and subtracts a linear drift from the head-on intervals surrounding each scan, constructs an ion beamâ€“beam deflection template, fits the remaining BPM motion, and exports per-step ROD corrections and diagnostic plots.

## Repository contents

| Path | Purpose |
| --- | --- |
| `residuals_fixed_gamma.ipynb` | End-to-end analysis, including the deflection calculation, fits, plots, and CSV exports. |
| `beam_beam_ion.py` | Standalone two-dimensional Gaussian field and ion-displacement functions; the notebook currently uses its **own inline copy**. See the sign-convention note below. |
| `set_scans_9212.csv` | Nominal scan steps, set positions, and time windows for fill 9212. |
| `9212.csv` | Time-stamped DOROS and arc BPM readings for fill 9212. |
| `linod_files/` | Example linear orbit-drift JSON outputs. |
| `resod_files_fill9212/` | Example per-step residual-correction CSV outputs. |
| `figures_fill9212/` | Example diagnostic figures (PNG and PDF). |

## Run the fill 9212 example

Use Python with Jupyter and the packages imported by the notebook:

```bash
python -m pip install jupyter numpy scipy matplotlib mplhep lmfit cycler
jupyter lab residuals_fixed_gamma.ipynb
```

In the notebook, set `PATH` in the **â€œread the file with the nominal valuesâ€** cell to the directory containing `set_scans_9212.csv` and `9212.csv`. It is a string used by concatenation, so include a trailing slash; for example, `PATH = './'` when Jupyter's working directory is the repository root. Then restart the kernel and run all cells in order. The notebook writes its outputs below `PATH`; rerunning it can overwrite existing figures and correction CSVs.

The input formats are:

- `set_scans_9212.csv`: **no header**; the notebook expects the columns, in order, `folder, scannr, subscan, axis, beam, typenr, step, timestart, timestop, reldisplacement, dispx, dispy, setX1, setX2, setY1, setY2`. Time values are Unix seconds and set positions are in mm. Rows with `subscan=Zero` define the surrounding head-on time windows.
- `9212.csv`: header followed by `Unix_Time` and left/right arc BPM and DOROS readings for X/Y and beams 1/2. BPM readings are in micrometres; the notebook converts them to mm internally.

## Analysis model

For each scan step, the notebook defines transverse separation as
`d = (setX1 - setX2, setY1 - setY2)` and adds the two beam covariance matrices to obtain the convolved profile. It evaluates a two-dimensional Gaussian Bassettiâ€“Erskine field `K(d, C1 + C2)`. For the same ion species in both beams, the kick magnitudes scale with
`2 N_opposing (ZÂ²/A) r_p / gamma_ion`: each beam uses the **opposing beam's ions per bunch**, calculated from the configured charge counts as `N_ions = N_charge / Z`. The notebook uses `gamma_ion = E_per_nucleon / m_p`.

With the notebook's separation and position convention, beam 1 uses `+K` and beam 2 uses `-K`. For each transverse plane `u`, it converts the kick angle `theta_u` to the arc/IP template and DOROS template with

```text
arc_u   = theta_u * beta_star / (2 tan(pi Q_u))
DOROS_u = theta_u * (beta_star / (2 tan(pi Q_u)) + doros_arm / 2)
```

All quantities in this calculation use mm for positions and lengths. The additional `mod` template simultaneously increases both bunch populations by 10% and transverse sizes by 10%; it is a **control variation**, not a measured systematic uncertainty.

The fit parameter `alpha` describes the nominal separation scale and `beta` the beamâ€“beam template amplitude. The notebook compares separate fits (`Sep`) with common length-scale (`CommonLS`), common beamâ€“beam (`CommonBB`), and combined common (`CommonLSBB`) fits, including subsets of vdM scans. Uncertainties for free `alpha` and `beta` parameters come from the `lmfit` covariance estimate. A fixed parameter has no fitted uncertainty; if a free parameter lacks one, the notebook stops and names the affected fit.

## Outputs and interpretation

- `linod_files/OrbitDrift_<fill>DOROS.json` and `OrbitDrift_<fill>arcBPM.json` contain the time windows and before/after-scan linear drift estimates. The checked-in example JSON filenames include an extra underscore after the fill number; the current notebook writes the names shown here.
- `resod_files_fill<fill>/` contains CSV corrections for each fit method and BPM source. Columns are `scan_name, plane, scan_point, correction_x, correction_y`; the exported corrections are in micrometres.
- `figures_fill<fill>/` contains BPM time traces, beamâ€“beam control curves, fitted length scales and amplitudes, and residual-versus-step plots.

The provided settings use a fixed covariance for Pb-208 in fill 9212, rather than measured covariance for each bunch crossing. Adapting the analysis to another AA system requires configuring its `Z`, `A`, beam populations, per-nucleon energy, transverse covariance, optics, and appropriately formatted scan/BPM inputs. The notebook does not separately model crossing-angle geometry in the deflection function.

**Sign-convention note:** `beam_beam_ion.py` currently assigns kick signs opposite to the inline `ion_bb_displacements` function in the notebook for the same `set(B1) - set(B2)` separation. Running the notebook uses the inline function. Check the coordinate convention against your BPM data before importing or substituting the standalone module.
