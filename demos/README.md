# Demonstrations of the Rayleigh integral and the angular spectrum method

This folder contains scripts that demonstrate how to use the Rayleigh integral and angular spectrum tools provided in this repository.

### Demos of the Rayleigh integral
- `demo_rayleigh.m`
  - Simple demonstration of the Rayleigh integral (one-dimensional cross-section of a pressure field).
- `demo_rayleigh_plane.m`
  - Demonstration of the Rayleigh integral (propagation to a plane)
- `demo_rayleigh_time.m`
  - This script compares a pressure field computed with the time-domain Rayleigh integral to a pressure field computed with the frequency-domain Rayleigh integral.

### Demos of the angular spectrum method
- `demo_angular_spectrum.m`
  - Basic demonstration of the angular spectrum method.
  - Corresponds to Fig. 6.1 of my thesis.
- `demo_find_angles.m`
  - This script demonstrates how to find the measurement plane orientation (Section 6.2.4 of my thesis).
  - Corresponds to Fig. 6.3 of my thesis.
  - Run `demo_rayleigh_plane.m` first to generate synthetic measurement data.
- `demo_angular_spectrum_volume.m`
  - This script demonstrates how the angular spectrum method can be used to rapidly generate the volume of the transducer output field by propagating to closely spaced planes.
- `demo_angular_spectrum_interpolation.m`
  - This script demonstrates that the angular spectrum method can also be applied to an inclined sensor plane through the use of interpolation.


### Comparison of propagation methods
- `comparison_rayleigh_angular_spectrum.m`
  - Comparison of a Rayleigh integral computation and an angular spectrum computation (one-dimensional cross section profile).
- `comparison_propagation_methods.m`
  - This script shows a comparison of propagation methods and source presentations as described in Section 6.2.3 of my thesis.
  - Compares Rayleigh integral, angular spectrum, and k-Wave.
  - This script requires installation of [k-Wave](https://github.com/PROTEUS-SIM/PROTEUS?tab=readme-ov-file#k-wave).
- `comparison_propagation_methods_plot.m`
  - Visualize and compare the data generated with `comparison_propagation_methods.m`.
  - Corresponds to Fig. 6.2 of my thesis.

### PROTEUS simulation settings for quick setup

- To enable a quick simulation setup for these demos, I have included a simulation setup file from PROTEUS in the current directory: `GUI_output_parameters.mat`.
- This settings file was created with the [graphical user interface of PROTEUS](https://github.com/PROTEUS-SIM/PROTEUS?tab=readme-ov-file#saving-settings-loading-settings-and-running-the-simulation).
- I have selected the default transducer.
- Note that this default transducer is a simplified model of the P4-1 transducer, based on manufacturer specifications, *not* the virtual transducer obtained with [main_pipeline.m](../main_pipeline.m).
- Settings modified in the PROTEUS GUI (*just for easy reference, no need to repeat, all the settings are stored in the MAT file*):
  1. `Transmit.NumberOfCycles = 1`
  2. `SimulationParameters.PointsPerWavelength = 4`
  3. `Geometry.Domain.Manual = true`
  4. `Geometry.Domain.Xmax = 0.010`
  5. `Medium.Tissue = 'Water'` (to load the properties of water)
  6. `Medium.Tissue = 'Custom'` (switched from `Water` to `Custom` to modify the following parameters)
  7. `Medium.AttenuationA = 0`
  8. `Medium.BonA = 0`
  9. `Medium.Inhomogeneity = 0`
 
### Note on new PROTEUS function (optional)

The branch `transducer-calibration` of the PROTEUS repository introduces the function [`update_sensor_fast`](https://github.com/PROTEUS-SIM/PROTEUS/blob/transducer-calibration/acoustic-module/update_sensor_fast.m), which is an upgrade of [`update_sensor`](https://github.com/PROTEUS-SIM/PROTEUS/blob/transducer-calibration/acoustic-module/update_sensor.m).
The new function is about two orders of magnitude faster, but this comes at the expense of a higher memory requirement.
- The maximum allowed RAM for the function is currently capped at 2 GiB, to prevent MATLAB from crashing on systems with limited RAM.
- If the memory estimate for the function exceeds the maximum value, `update_sensor_fast.m` redirects to `update_sensor.m`. This will make these scripts much slower:
  - [`comparison_propagation_methods.m`](comparison_propagation_methods.m)
  - [`demo_angular_spectrum_interpolation.m`](demo_angular_spectrum_interpolation.m)
- If sufficient RAM is available on your system, increase the value of `maxMemory` in `update_sensor_fast` to at least 4 GiB.

