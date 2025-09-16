# Figure formatting and exporting

This folder contains scripts for recreating the figures of Chapter 6 in my thesis.

These scripts only perform figure formatting and exporting.

The base plots are generated with:
- Fig. 1:  `demos/demo_angular_spectrum.m`
- Fig. 2:  `demos/comparison_propagation_methods_plot.m`
- Fig. 3:  `demos/demo_find_angles.m`
- Fig. 4:  no numerical data
- Fig. 5:  `main_pipeline.m > rayleigh_propagation.m`
- Fig. 6:  `main_pipeline.m > estimate_model_parameters.m` and `compute_transmit_impulse_response.m`
- Fig. 7:  `main_pipeline.m > determine_virtual_receiver_orientation.m` and `compute_receive_impulse_response.m`
- Fig. 8:  `analysis/validation.m`
- Fig. 9:  `analysis/analysis_modes.m`
- Fig. 10: `analysis/analysis_phase.m`

The base plots are stored in the folder `fig`.
Exported figures are stored in `exports`.

`CustomColorMap.mat` is a copy from [Super-Resolved Microbubble Localization in Single-Channel Ultrasound RF Signals Using Deep Learning](https://github.com/MIAGroupUT/SRML-1D/tree/main/RF_simulator/ColorMapMaker).
