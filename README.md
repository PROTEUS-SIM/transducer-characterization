# Characterization and virtualization of a medical ultrasound transducer

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](LICENSE)

Nathan Blanken

This repository contains the code for the preprint
[arXiv:2509.22090](https://arxiv.org/abs/2509.22090)
(doi: [10.48550/arXiv.2509.22090](https://doi.org/10.48550/arXiv.2509.22090)).

## System requirements

- Installation of MATLAB is required. The use of MATLAB 2019b or later is recommended.
- Although the characterization pipeline can be run on a CPU, the use of a GPU with at least 8 GB of dedicated memory is recommended to speed up the Rayleigh integral computations.
- To run the code on GPU, the Parallel Computing Toolbox™ from MATLAB is required.


## Installation

### Repository

Clone the current repository with git.
```
git clone https://github.com/NathanBlanken/transducer-characterization.git
```

### PROTEUS acoustic module

This repository depends on the latest functionality from the acoustic module of the PROTEUS toolbox.
This functionality currently exists in the branch
[`transducer-calibration`](https://github.com/PROTEUS-SIM/PROTEUS/tree/transducer-calibration) of the PROTEUS repository.
This branched will be merged with the
[main](https://github.com/PROTEUS-SIM/PROTEUS) branch in the next version of PROTEUS.

1. Clone PROTEUS with git:
   ```
   git clone https://github.com/PROTEUS-SIM/PROTEUS.git
   ```

3. On your computer, switch to the branch `transducer-calibration`:
   ```
   git checkout transducer-calibration
   ```

**Note:** a [full](https://github.com/PROTEUS-SIM/PROTEUS?tab=readme-ov-file#installation-of-the-simulator)
installation of PROTEUS, including k-Wave, geometry data, flow solver, and vtkToolbox, is **not** required
if you are only going use it for transducer characterization.

**Note (optional):** installation of [k-Wave](https://github.com/PROTEUS-SIM/PROTEUS?tab=readme-ov-file#k-wave) is required
if you intend to reproduce Fig. 2 or Fig. 8 of the arXiv preprint.

**Note (optional):** the branch `transducer-calibration` of the PROTEUS repository introduces the function [`update_sensor_fast`](https://github.com/PROTEUS-SIM/PROTEUS/blob/transducer-calibration/acoustic-module/update_sensor_fast.m), which is an upgrade of [`update_sensor`](https://github.com/PROTEUS-SIM/PROTEUS/blob/transducer-calibration/acoustic-module/update_sensor.m). To unlock the full potential of this new function, increase the value of `maxMemory` to a higher value that is still supported by your system's RAM. See [`demos`](demos) for further details.

### Path variables

After installation of PROTEUS, navigate back to the current repository ([`transducer-characterization`](.)).

Modify the variable `PATHS.AcousticModule` in the file [`path_setup_characterization.m`](path_setup_characterization.m)
to match the location of the PROTEUS acoustic module on your computer.
```
PATHS.AcousticModule = '/home/user/PROTEUS/acoustic-module'; % Replace by your own path
```

To add all required directories to the MATLAB path, run the file from the MATLAB command window:
```
path_setup_characterization
```
To remove those directories again, run:
```
path_setup_characterization('rmpath')
```


### Experimental data

Download the experimental data of the P4-1 transducer from
[Zenodo](https://doi.org/10.5281/zenodo.17095584)
to the folder `data`.

With this experimental data, you will be able to reproduce all results in the arXiv preprint.
Reproducing these results is recommended as a tutorial.

**Note:** the folder `PROTEUS-I` only needs to be downloaded,
if you want to reproduce Fig. 8 of the arXiv preprint.
This folder contains the data from Fig. 4 of [PROTEUS Part I](https://ieeexplore.ieee.org/document/10597664).

## Overview of the repository

- 📂 [`analysis`](analysis) Analysis of the characterization results.
  Corresponds to Sections III C, III D, and IV of the arXiv preprint.
- 📂 `angular-spectrum` Tools related to the angular spectrum method.
- 📂 `data` Experimental data.
- 📂 [`demos`](demos) Demonstrations of the functions in 📂 `angular-spectrum`, and 📂 `rayleigh-integral`.
  These demonstrations correspond to Section II of the arXiv preprint.
- 📂 [`figures`](figures) Figure formatting scripts used for the figures in the arXiv preprint.
  The readme in this folder contains an overview of how each figure was generated.
- 📂 `rayleigh-integral` Tools related to the Rayleigh integral.
- 📂 `results` Output folder for characterization results
    (impulse responses and transducer model parameters)
- 📂 `utilities` Helper functions.


### Demonstration of the full characterization pipeline

The script [`main_pipeline.m`](main_pipeline.m) in the root directory is
a demonstration of the full characterization pipeline.
It is an implementation of Fig. 4 in the arXiv preprint.
The other functions in the root folder are called by `main_pipeline.m` and
correspond to one or multiple steps in Fig. 4.

The file `main_pipeline.m` returns the transmit and receive impulse responses
and the transducer properties in a format that can be directly loaded into PROTEUS.
See [TransducerGUI](https://github.com/PROTEUS-SIM/PROTEUS/blob/main/documentation/TransducerGUI.md).

For the an interactive experience of the demonstration in `main_pipeline.m`,
it is recommended to use the checkpoints in the file,
which allow you to run one part of the pipeline at a time.
Run `help main_pipeline.m` for further details.

### Other demonstrations

If you are not so familiar with the theory presented in Section II of the arXiv preprint,
you are encouraged to check out the demos in [`demos`](demos) first.

### Coordinate system

The tools in this repository use a Cartesian coordinate system with the *z*-axis pointing along the axis of propagation.
By contrast, PROTEUS currently uses a Cartesian coordinate system with *x*-axis pointing along the axis of propagation.
The folder `utilities` contains tools for switching smoothly between the two systems.
See [`demos`](demos) for further information.

### Applying the pipeline to new transducer characterization data

Note that `main_pipeline.m` is a prototype framework, not an out-of-the-box tool
for transducer characterization.
Applying `main_pipeline.m` to new data requires fine-tuning of the computational steps,
in particular the parameters at the top of each of these files:
 - `compute_receive_impulse_response.m`
 - `determine_scan_plane_orientation.m`
 - `determine_virtual_receiver_orientation.m`
 - `estimate_model_parameters.m`

The pipeline was demonstrated on a transducer with a flat, rectangular aperture.
Applying the pipeline the pipeline to an arbitrary transducer geometry will require more substantial modification of the code.

For example, to find the radius of curvature for a curved transducer, `estimate_model_parameters.m` has to be modified.
Subsequently, the file `rayleigh_propagation` has to be modified to set up curved source/sensor structures.

For transducers with a non-rectangular aperture, `determine_scan_plane_orientation.m` has to be modified.
(This involves merely a simplification in the case of a circular aperture.)

Nonetheless, the tools in 📂 `rayleigh-integral` are already fully compatible with an arbitrary geometry.

### Preparing hydrophone data in the correct format

Hydrophone scan data must be stored as a three-dimensional array with dimensions Nx-by-Ny-by-Nsamples,
where Nx-by-Ny is the size of the size of the scan grid and the third dimension corresponds to the time-domain hydrophone voltage data.

Additionally,
a MATLAB structure `medium`, with fields `sound_speed` and `density`, and
a MATLABstructure `Grid` in [PROTEUS format](https://github.com/PROTEUS-SIM/PROTEUS/blob/main/acoustic-module/define_grid.m) must be provided.

Compare your data format to the data format in `Transmit.mat` from [Zenodo](https://doi.org/10.5281/zenodo.17095584).

**Note for collaborators at the University of Twente:**
Collaborators acquiring data with the [calibration setup](https://github.com/NathanBlanken/calibration-setup) at the Physics of Fluids group
can use the script [`reorganize_data`](utilities/reorganize_data.m)
to convert the experimental data to the correct format.


## License and citation

Copyright (C) 2025 Nathan Blanken

This code repository is licensed under the GNU Lesser General Public License v3.0 (LGPL-3.0).
See the [LICENSE](./LICENSE) file for further details.
GitHub will lists an overview of permissions, limitations, and conditions on top of the page.

If you use (parts of) the code, please consider citing the arXiv preprint:

N. Blanken, M. Versluis, and G. Lajoinie,
_Characterization and virtualization of a medical ultrasound transducer_ (2025),
arXiv:2509.22090.

```
@misc{Blanken2025,
    author = {Nathan Blanken and Michel Versluis and Guillaume Lajoinie},
    doi    = {10.48550/arXiv.2509.22090},
    eprint = {2509.22090},
    title  = {Characterization and virtualization of a medical ultrasound transducer},
    url    = {https://arxiv.org/abs/2509.22090},
    year   = 2025,
    archivePrefix = {arXiv},
    primaryClass  = {physics.med-ph},
}

```
