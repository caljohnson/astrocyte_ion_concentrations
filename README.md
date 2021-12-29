# Astrocyte Ion Concentrations
Astrocyte model to examine the level of astrocyte calcium activity needed to have a significant impact on extracellular concentrations, and hence, the excitability of nearby neurons.

## Folders

* ./ - Main codes to run simulations, parameter sweeps, generate data and figures
* ./EIF_model - Codes to compare EIF neuron and WB neuron models
* ./WB_neuron_xpp_ode - XPPAUT scripts for generating bifurcation diagrams for the WB neuron model
* ./src - Source codes and functions that main codes run off of
* ./testing_existing_models - Codes to test and fit individual channel/transporter models to relevant data

<details>
<summary>Model Fitting & Validation Figures</summary>

* Figure 2 (Astrocyte model and experimental \cite{kirischuk2012sodium} (a) internal calcium response and (b) internal sodium response to a 1 mM bath-application of (external) glutamate)
    * Run "Glubath_compare_with_Data.m"
