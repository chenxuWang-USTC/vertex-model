# Vertex Model for Organoid Mechanics

This repository contains MATLAB models for simulating the mechanical behavior, growth, and morphological evolution of organoids.

The repository includes:

- a two-dimensional vertex/trapezoidal model for organoid cross-sections;
- a simplified energy model for rapidly exploring organoid growth trends.

The models investigate how cell-surface tension, internal pressure, viscous dynamics, extracellular-matrix (ECM) interactions, ablation, and osmotic perturbations influence organoid morphology.

> Reference repository: [chenxuWang-USTC/vertex-model](https://github.com/chenxuWang-USTC/vertex-model)

## Software and demo-data availability

The repository provides MATLAB source code rather than a compiled standalone executable. No external dataset is required for the supplied demonstration: the examples generate simulated states from the initial conditions and parameters included in the code. Users who want to model experimental samples can replace those initial conditions and parameters as described in [Using the software with your own data](#using-the-software-with-your-own-data).

## Repository structure

```text
vertex-model/
├── 2D vertex model/
│   ├── Morphological Evolution of Organoids/
│   ├── ablation_healing/
│   ├── osm_change/
│   └── README.md
├── brief energy model/
│   ├── change_alpha.m
│   └── README.md
└── LICENSE
```

### `2D vertex model`

The two-dimensional model represents a cross-section of an organoid using vertex/trapezoidal approximations and viscous dynamics. It contains three groups of simulations:

1. global morphological evolution;
2. laser ablation and subsequent healing;
3. response to osmotic-pressure changes.

### `brief energy model`

The brief energy model provides a faster and simpler estimate of organoid growth. The `change_alpha.m` script illustrates the relationship between cell number and apical actin density and generates the video `cell_grow.avi`.

## Requirements and tested versions

- **Required software:** MATLAB R2018b or later.
- **Tested environment:** MATLAB R2023b on 64-bit Windows 11.
- **Other operating systems:** The source code is expected to run on recent 64-bit versions of Windows, macOS, and Linux supported by MATLAB, but these combinations have not all been individually verified.
- **MATLAB dependencies:** No third-party MATLAB package or non-standard toolbox is listed as a required dependency.
- **Hardware:** No non-standard hardware, GPU, or laboratory instrument is required.
- **Recommended desktop configuration:** A 4-core CPU, 8 GB RAM, and at least 500 MB of free disk space for the source code and generated output.

The current repository does not include an automated test suite. MATLAB R2018b is the minimum supported version; the environment listed above is the version used for validation. Small numerical or graphical differences may occur between MATLAB releases, operating systems, hardware, and parameter settings.

A run is considered successful when:

- the selected script completes without an error;
- the state variables do not contain unexpected `NaN` or `Inf` values;
- the expected geometry figures, time-series plots, or video are generated;
- the qualitative response agrees with the perturbation applied to the model.

## Installation

**Typical installation time:** approximately **10 minutes** on a normal desktop computer. This estimate includes downloading the repository, opening MATLAB, adding the folders to the MATLAB path, and checking the helper functions. MATLAB installation time is not included.

### 1. Clone or download the repository

```bash
git clone https://github.com/chenxuWang-USTC/vertex-model.git
cd vertex-model
```

Alternatively, use **Code → Download ZIP** on GitHub and extract the archive.

### 2. Start MATLAB in the repository root

In MATLAB, select the cloned `vertex-model` directory as the current folder. Then add all model directories to the MATLAB path:

```matlab
repoRoot = pwd;
addpath(genpath(repoRoot));
```

To keep the path after MATLAB restarts, optionally run:

```matlab
savepath;
```

If `savepath` is not permitted, omit it and rerun `addpath(genpath(repoRoot))` in each new MATLAB session.

### 3. Verify the helper functions

The 2D vertex model depends on helper functions such as:

- `intital.m`;
- `quad_area.m`;
- `compute_total_energy.m`.

The spelling `intital.m` is retained from the source repository.

Check that MATLAB can locate these files:

```matlab
which intital
which quad_area
which compute_total_energy
```

Each command should return a valid file path. If a path is missing, make sure the corresponding `.m` file is present and that the repository has been added recursively to the MATLAB path.

## Demo

The fastest demonstration uses the brief energy model. It does not require an external experimental dataset: the script uses the parameters and simulated initial conditions supplied with the source code.

**Expected demo runtime:** approximately **10 minutes or less** on a normal desktop computer. The precise runtime depends on the MATLAB release, CPU, parameter values, iteration count, plotting speed, and video-writing settings.

### Demo instructions

1. Start MATLAB and select the repository root as the current folder.
2. Add the complete repository to the MATLAB path:

   ```matlab
   repoRoot = pwd;
   addpath(genpath(repoRoot));
   ```

3. Run the brief energy-model demo:

   ```matlab
   cd(fullfile(repoRoot, 'brief energy model'));
   change_alpha;
   ```

4. Confirm that the script finishes without an error and that `cell_grow.avi` has been created:

   ```matlab
   assert(isfile('cell_grow.avi'), ...
       'Demo failed: cell_grow.avi was not created.');
   implay('cell_grow.avi');
   ```

### Expected demo output

- MATLAB displays the simulated organoid-growth process.
- A video named `cell_grow.avi` is written to `brief energy model` or to the active output directory.
- The video shows the modeled relationship between cell number and apical actin density.
- A successful demonstration completes without an undefined-function, missing-file, `NaN`, or `Inf` error.

The 2D vertex-model modules described below provide additional, more detailed demonstrations. Their runtime may exceed 10 minutes if the mesh, number of iterations, output frequency, or parameter scan is increased.

## How to run the models

The exact top-level script name can differ among the subdirectories. In each module, run the simulation script rather than an individual helper function. Before a full simulation, inspect the parameter block at the beginning of the selected script.

### 1. Morphological evolution of organoids

This module studies the global equilibrium shape of an organoid under different mechanical parameters.

Change to the module directory:

```matlab
cd(fullfile(repoRoot, '2D vertex model', ...
    'Morphological Evolution of Organoids'));
```

Open and run the top-level simulation script in this directory. Relevant parameters may include:

- apical, basal, and lateral surface tensions;
- internal pressure;
- cell or tissue viscosity;
- ECM stiffness or coupling;
- integration time step;
- total simulation time and convergence criteria.

#### Expected result

- MATLAB displays the organoid cross-section or cell geometry.
- The initial geometry evolves toward a mechanically balanced shape.
- Depending on the script, MATLAB may also display total energy, area, force, displacement, or convergence curves.
- Changing pressure or surface tension should alter the organoid size, curvature, or degree of deformation.
- Near equilibrium, the displacement or energy change between successive steps should decrease.

Exact values depend on the initial geometry and chosen parameters; the expected result is therefore qualitative rather than a fixed numerical value.

### 2. Ablation and healing

This module simulates laser ablation of a single cell or a group of cells, followed by mechanical relaxation and topological healing.

Change to the module directory:

```matlab
cd(fullfile(repoRoot, '2D vertex model', 'ablation_healing'));
```

Open the top-level simulation script and check the ablation target, ablation time, and post-ablation simulation duration before running it.

#### Expected result

- The selected cell, junction, or local mechanical contribution is removed or weakened.
- Neighboring vertices show an immediate recoil or redistribution of force.
- Surrounding cells rearrange during the healing stage.
- The ablated region may partially or completely close, depending on the selected parameters.
- Figures show the geometry before ablation, immediately after ablation, and during or after healing.

### 3. Osmotic-pressure change

This module explores the dynamic response of an organoid to a change in internal hydrostatic pressure, representing an osmotic shock.

Change to the module directory:

```matlab
cd(fullfile(repoRoot, '2D vertex model', 'osm_change'));
```

Before running the top-level script, check the initial pressure, changed pressure, perturbation time, integration time step, and simulation duration.

#### Expected result

- Increasing the internal pressure generally produces organoid expansion.
- Decreasing the internal pressure generally produces contraction.
- The geometry relaxes toward a new mechanical state after the pressure change.
- The generated figures or time-series data may show changes in cross-sectional area, shape, displacement, pressure, or total energy.

The precise response depends on the sign and magnitude of the pressure change and on the mechanical parameters of the organoid and ECM.

### 4. Brief energy model

This model gives a rapid, approximate picture of organoid growth and the relationship between cell number and apical actin density.

Run:

```matlab
cd(fullfile(repoRoot, 'brief energy model'));
change_alpha;
```

#### Expected result

- MATLAB displays the modeled growth process.
- The script creates `cell_grow.avi` in the current working directory.
- The video illustrates the relationship between cell number and apical actin density.

Play the generated video in MATLAB with:

```matlab
implay('cell_grow.avi');
```

If `cell_grow.avi` already exists, rename or move it before running the script if the previous output must be preserved.

## Using the software with your own data

These programs are primarily forward simulations rather than a fixed experimental-data analysis pipeline. Consequently, “your own data” normally means supplying your own initial geometry and experimentally estimated mechanical parameters.

### 1. Prepare the inputs

Convert experimental measurements to the variables and units expected by the selected script. Depending on the module, relevant inputs can include:

- vertex coordinates or the initial organoid cross-section;
- cell connectivity or cell-neighbor information;
- apical, basal, and lateral tensions;
- internal pressure or a pressure-change schedule;
- tissue viscosity and ECM stiffness/coupling;
- ablated cell or junction indices;
- time step, total simulation time, and output interval.

The repository does not define a universal CSV or MAT-file schema for every module. Match the dimensions, indexing convention, and units used by the example variables in the corresponding script and initialization function.

### 2. Create a separate parameter copy

Keep the original example unchanged. Copy the top-level script or its parameter block, give the copy a descriptive name, and change only the input geometry or parameter values. This makes the supplied demonstration reproducible and prevents accidental loss of the defaults.

### 3. Update the initialization

If custom geometry is required, edit a copy of `intital.m` or replace the initial-coordinate and connectivity variables after initialization. Confirm that:

- vertex indices start at 1, as required by MATLAB;
- all referenced vertices and cells exist;
- coordinate and parameter units are mutually consistent;
- polygon ordering is consistent and does not create self-intersections;
- areas and lengths are positive before starting the simulation.

### 4. Run and validate

Run the selected top-level script from its own module directory. First use a short simulation or reduced iteration count. Inspect the initial geometry, check the MATLAB Command Window for warnings, and verify that the state variables remain finite before starting a full parameter scan.

Because input formats differ between modules, document any conversion script and include the converted input file, units, and parameter values with the resulting figures or videos.

## Quick verification

After installation, use the following checklist:

1. Start MATLAB in the repository root.
2. Run `addpath(genpath(pwd))`.
3. Confirm that `which intital`, `which quad_area`, and `which compute_total_energy` return valid paths.
4. Run one top-level script from the desired 2D vertex-model module.
5. Confirm that figures are generated and that no undefined-function or missing-file error occurs.
6. Run `change_alpha.m` from `brief energy model`.
7. Confirm that `cell_grow.avi` is created and can be played.

## Troubleshooting

### MATLAB reports `Undefined function or variable`

Re-add the repository recursively to the MATLAB path:

```matlab
cd(repoRoot);
addpath(genpath(repoRoot));
rehash;
```

Check that the spelling and capitalization in the error message match the file name on disk.

### MATLAB cannot find a data or helper file

Some scripts may use relative paths. Run the script from its own directory, or update the path in the script to use `fullfile`. Do not move individual `.m` files away from the module unless all relative dependencies are updated.

### The simulation returns `NaN` or `Inf`

- Reduce the integration time step.
- Restore the default parameter values.
- Check that pressure, viscosity, stiffness, tension, and area parameters use physically valid signs and units.
- Modify one parameter at a time to identify the unstable combination.

### The geometry does not change

Check that:

- the total simulation time and number of iterations are nonzero;
- the driving pressure or tension difference is nonzero;
- the viscosity or damping is not so large that motion is visually negligible;
- the plotted frame is updated inside the simulation loop.

### `cell_grow.avi` is not created

- Use `pwd` to confirm the output directory.
- Confirm that MATLAB has write permission in that directory.
- Close any program that already has `cell_grow.avi` open.
- Check whether an existing file is being overwritten or locked.

## Reproducibility

If a model uses random initialization, set a random seed before running it:

```matlab
rng(1, 'twister');
```

For reproducible simulations, record:

- the Git commit or release;
- the MATLAB release and operating system;
- all modified model parameters;
- the random seed, if applicable;
- the directory and script used as the simulation entry point.

### Reproducing a saved result

To reproduce a quantitative result or manuscript figure:

1. Check out the same Git commit or archived release.
2. Use the recorded MATLAB release and operating system where possible.
3. Restore the original input geometry and complete parameter set.
4. Set the recorded random seed before initialization.
5. Run the same top-level script from the same module directory.
6. Save the MATLAB Command Window output, generated figures/data, and the parameter file together.
7. Compare both the qualitative morphology and the reported numerical quantities; do not rely only on the appearance of a plot.

## Software-submission checklist coverage

This README provides the items normally requested for a software/code submission:

| Required item | Location in this README |
| --- | --- |
| Source code availability | Repository structure and installation sections |
| Software and operating-system requirements | Requirements and tested versions |
| Versions tested | Requirements and tested versions |
| Non-standard hardware requirements | Requirements and tested versions |
| Installation instructions | Installation |
| Typical installation time | Installation |
| Demo instructions | Demo |
| Expected demo output | Demo |
| Expected demo runtime | Demo |
| Instructions for use with custom data | Using the software with your own data |
| Reproduction guidance | Reproducibility |

## Scope and limitations

- The 2D vertex model represents an organoid cross-section rather than a full three-dimensional organoid.
- The brief energy model is intended for rapid qualitative exploration and does not replace the more detailed vertex-model dynamics.
- Model output depends on parameter choices and should be interpreted as a computational prediction rather than direct experimental evidence.

## License

This repository is distributed under the [MIT License](LICENSE).

Copyright © 2026 organoid.

## Citation

No formal citation information is currently provided in the reference repository. If the code contributes to a publication, cite the repository URL and the exact commit or archived release used to generate the results.
