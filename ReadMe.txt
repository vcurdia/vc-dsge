-*- mode: org -*-

*vcDSGE toolbox*

v1.0.0

vcDSGE is a toolbox to simulate and estimate DSGE models.

Created: January 21, 2016 by Vasco Curdia

Copyright 2016-2017 Vasco Curdia


* License
vcTools is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

vcTools is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
vcTools. If not, see <http://www.gnu.org/licenses/>.


* Requirements

** Matlab
Matlab R1016b with following toolboxes
- Symbolic Toolbox
- Statistical Toolbox
- Optimization Toolbox

** LaTeX
LaTeX is used to compile certain documents.

** vcTools by Vasco Curdia
Miscellaneous scripts, functions and objects.
Provided with vcDSGE.

** gensys by Chris Sims
DSGE rational expectations equilibrium solver.
Available at http://sims.princeton.edu/yftp/gensys/

** optimize by Chris Sims
Set of optimization and non-linear equation solution scripts and functions.
Available at http://dge.repec.org/codes/sims/optimize/

** KF by Chris Sims
Kalman filter and smoother functions.
Available at http://sims.princeton.edu/yftp/Times09/KFmatlab/


* Description

The vcDSGE toolbox is organized around the following key objects and one
example.

** DSGE.Model
Contains all model components, including variables, parameters, equations.
It also sets the methods to:
- solve and evaluate the DSGE model
- simulate the model, including:
  - IRF: impulse response functions
  - VD: variance decomposition
  - States: simulate the DSGE states over sample period
  - SD: shock decomposition of states through sample

It does not have to be and estimated model. It can be a simple parameterized
model.

For more details type the following inside Matlab:
doc DSGE.Model

** DSGE.Prior
Object with all components regarding the DSGE parameters prior. It includes
methods to:
- report its properties
- draw from the prior
- evaluate pdf/log-pdf

For more details type the following inside Matlab:
doc DSGE.Prior

** DSGE.Posterior
Object with DSGE parameter posterior components. Includes methods to:
- Evaluate posterior pdf/log-pdf
- maximize posterior
- generate MCMC sample
- analyze MCMC sample convergence
- report parameters from MCMC sample

For more details type the following inside Matlab:
doc DSGE.Posterior

** DSGE.Data
Object with data and data properties to use throughout the DSGE estimation and
simulation.

For more details type the following inside Matlab:
doc DSGE.Data

** Example
In the subfolder Example there are scripts and data to show how to use the main
features of vcDSGE:

*** setupDSGE
setup a simple DSGE model, including its prior, posterior and data. It also
performs some simulations using the prior distribution.

*** estimateDSGE
estimates the DSGE model setup in the previous code and performs simulations
using posterior sample.
