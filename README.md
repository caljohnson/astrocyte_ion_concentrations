# Astrocyte Ion Concentrations
Astrocyte model to examine the level of astrocyte calcium activity needed to have a significant impact on extracellular concentrations, and hence, the excitability of nearby neurons.

## Folders

* ./ - Main codes to run simulations, parameter sweeps, generate data and figures
* ./EIF_model - Codes to compare EIF neuron and WB neuron models
* ./WB_neuron_xpp_ode - XPPAUT scripts for generating bifurcation diagrams for the WB neuron model
* ./src - Source codes and functions that main codes run off of
* ./testing_existing_models - Codes to test and fit individual channel/transporter models to relevant data

<details>
<summary>How to Generate Figures</summary>

* Figure 2 (Astrocyte model and experimental (a) internal calcium response and (b) internal sodium response to a 1 mM bath-application of (external) glutamate)
    * Run "Glubath_compare_with_Data.m"
* Figure 3 (Astrocyte model (a) external potassium clearance and (b) sodium efflux in response to elevated external potassium conditions (9 mM))
    * Run "astrocyte_modulatory_role_elevatedK_tuning.m"
* Figure 4 (Bifurcation diagrams for the WB neuron model)
    * Run "wb_neuron.ode" in XPP.  
    *To generate the bifurcation diagram, start at a high Iapp, run to steady state in XPP, and launch AUTO.
        *Use Iapp as the bifurcation parameter and run with negative parameter steps to find the Hopf bifurcation point.
    *To generate the frequency/period plots, generate the bifurcation diagram then change the plot type to frequency/period.
* Figures 5-6 (Simulated ENa,EK stability diagrams for the WB neuron)
    * Run "remote_WB_neuron_ena_ek_loop.m" to generate the data
    * Run "figs_remote_WB_neuron_ena_ek_loop.m" to generate the figures from the data
* Figure 7 (3D stability diagram)
    * Run "remote_WB_neuron_ena_ek_loop.m" to generate the data
    * Run "figs_remote_WB_neuron_ena_ek_loop_3dbifdiag.m" to generate the figures from the data
* Figure 8 (Astrocyte model simulated with and without the IP3-mediated calcium transient in response to a glutamate stimulus)
    * Run "test_full_system_with_Ca_forcing_catransient_andGlu.m"
* Figure 9 and 11 (Astrocyte trajectory overlaid on stability diagrams)
    * Run "test_full_system_with_Ca_forcing_catransient_andGlu.m" to generate astrocyte-effect data
    * Run "test_full_system_with_Ca_forcing_elevatedK_catransient" to generate astrocyte-effect data
    * Run "remote_WB_neuron_ena_ek_loop.m" to generate the WB neuron model data
    * Run "figs_remote_WB_neuron_3dbifdiag_withastrocyte_trajectories.m" to generate the figures
* Figure 10 (Astrocyte model simulated with and without the IP3-mediated calcium transient in response to elevated external potassium [K+]e(0) = 9 mM)
    * Run "test_full_system_with_Ca_forcing_elevatedK_catransient"

