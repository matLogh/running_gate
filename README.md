# running_gate

RunningGate is a small C++ library (using the ROOT framework) implementing the "running gate" technique
for 2D histograms (see: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.58.771). It creates sequential
1D projections from a `TH2` matrix, fits peaks in each projection and produces evolution graphs (multigraphs)
for peak areas vs gate center energy.

Main capabilities
- Create projections from a `TH2` and evaluate peak areas per gate
- Multiple fitting backends: event-counting, Gaussian fitter; integrates the bundled `theuerkauf_fitter` module
- Save projections, fits and canvases into a ROOT file (fits are written so opening the file in a new session recreates plots)
- Helpers: convert `TGraphErrors` → `TH1F`, find X at maximum Y, produce multigraph with confidence intervals

Requirements
- ROOT (tested with ROOT 6.24)
- C++17 toolchain

Build
1. Configure and build with CMake (recommended):

```bash
mkdir build; cd build
cmake ../
cmake -j4
```

2. The build links example data to the build directory to run the provided example executable **running_gate**.

Basic usage (in code)

```cpp
// create RunningGate from a TH2* (your 2D matrix)
RunningGate rg(matrix);

// register peaks to follow (choose fitter)
rg.SetCoincidencePeak_gausFit(peak_center1, peak_start1, peak_end1, bcg_left_start1, bcg_left_end1, bcg_right_start1, bcg_right_end1, background_polynomial_order1);
rg.SetCoincidencePeak_countEvents(peak_center1, peak_start1, peak_end1, bcg_left_start1, bcg_left_end1, bcg_right_start1, bcg_right_end1, background_polynomial_order2);
// recommended but optional save output
rg.SetSaveFileName("running_gate_output.root");

// run analysis
rg.RunGate(energy_gate_start, energy_gate_end, bin_gate_width);
```

Notes & tips
- It is recommended to use save file option for diagnostics. The file contain all gate projections, count/fit of specific coincidence peaks for each gate projection, all graphs, canvases and *available only in file* is the analysis canvas of the gaussian fit (Theuerkauf peak model, gaussian with tails) of the resulting running gate spectrum
- If running headless (no GUI) use `gROOT->SetBatch(kTRUE)` to suppress windows.


Examples
- `ExampleCode.cpp` shows a runnable example and how to call `running_gate`.


License
See repository files for license information.
