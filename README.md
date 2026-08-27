# vc-dsge toolbox

An object-oriented MATLAB toolbox to simulate and estimate Dynamic Stochastic
General Equilibrium (DSGE) models, including Bayesian estimation via MCMC.

Created January 21, 2016 by Vasco Cúrdia.

## Requirements

### MATLAB

Developed and tested on MATLAB R2024a. Other releases are untested and may
work. The following toolboxes are required:

- Symbolic Math Toolbox
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

### Other dependencies

| Dependency | Purpose | Source |
| --- | --- | --- |
| vc-tools (Vasco Cúrdia) | Miscellaneous scripts, functions and objects | [github.com/vcurdia/vc-tools](https://github.com/vcurdia/vc-tools) |
| gensys (Chris Sims) | DSGE rational expectations equilibrium solver | <http://sims.princeton.edu/yftp/gensys/> |
| optimize (Chris Sims) | Optimization and non-linear equation solvers | <http://dge.repec.org/codes/sims/optimize/> |
| KF (Chris Sims) | Kalman filter and smoother | <http://sims.princeton.edu/yftp/Times09/KFmatlab/> |

LaTeX is used to compile certain output documents, and the example sets the
default MATLAB text interpreter to `latex`.

## Setup

Add this toolbox and its dependencies to the MATLAB path. `example/setpath.m`
shows one such arrangement, expecting the following layout under `~/matlab`:

```
~/matlab/
  vc-dsge/          <- this toolbox
  vc-tools/
  sims/
    gensys/
    kf/
    optimize/
```

Adjust `pathBase` and `pathList` in `example/setpath.m` to match your own layout.

## Quick start

The `example/` directory contains a complete worked model:

| File | Purpose |
| --- | --- |
| `setupdsge.m` | Sets up a simple DSGE model with its prior, posterior and data, then runs simulations from the prior |
| `estimatedsge.m` | Estimates that model and runs simulations from the posterior sample |
| `setpath.m` | Adds the toolbox and dependencies to the MATLAB path |
| `Data_1987q3_2009q3.csv` | Sample data |

Run `setupdsge` first, then `estimatedsge`.

## Description

The toolbox is organized around classes in the `DSGE` package namespace. All
derive from `matlab.mixin.Copyable`, so instances have value-like copy
semantics via `copy()`.

### `DSGE.Model`

The central class. Holds all model components — variables, parameters,
equations, prior, posterior and data — and provides the methods that do the
work:

- **Solve and evaluate** — `solveree`, `genmats`, `expandparam`
- **Simulate** — `irf` (impulse responses), `vd` (variance decomposition),
  `states` (states over the sample), `sd` (shock decomposition), `sim`
- **Prior** — `initializeprior`, `priordraw`, `priorpdf`, `priorlpdf`,
  `makepriorsample`, `analyzeprior`
- **Posterior** — `initializepost`, `postlpdf`, `postdraw`, `maxpost`,
  `maxpostreport`, `maxpostchoosebest`, `pickmaxpost`, `analyzepost`
- **MCMC** — `mcmc`, `mcmcchain`, `mcmcconvergence`, `mcmcredux`,
  `calibratemcmc`, `loadmcmcdraws`
- **Output** — `figpanels`

A model need not be estimated; it can be a simple parameterized model used only
for simulation.

### `DSGE.MCMC`

Encapsulates an MCMC run over the model posterior, including its configuration
and the `run` method that generates the sample.

### `DSGE.Data`

Data and data properties used throughout estimation and simulation.

### `DSGE.Var` and `DSGE.Param`

Model variables and parameters, including the parameter prior specification.

### `DSGE.Options.Table`

Options controlling formatted table output.

For details on any class or method, use MATLAB's help, e.g.:

```matlab
doc DSGE.Model
help DSGE.Model.irf
```

## License

vc-dsge is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

vc-dsge is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
vc-dsge. If not, see <http://www.gnu.org/licenses/>.

Copyright 2016-2026 Vasco Cúrdia
