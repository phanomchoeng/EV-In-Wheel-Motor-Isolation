# EV In-Wheel Motor Isolation

MATLAB source code accompanying the manuscript:

**“Mitigating Unsprung Mass in an Electric Vehicle through Passive and Active Isolation of a Wheel-Integrated Motor”**

## Description

This repository contains the MATLAB source code used to reproduce the numerical analyses presented in the manuscript.

The study compares three mass-matched electric-vehicle suspension configurations:

1. Centralized-drive electric vehicle
2. Rigid wheel-mounted motor
3. Passively and actively isolated wheel-integrated motor

The numerical analyses include:

- Mass-matched suspension models
- Passive motor-isolation parameter search
- Motor-mass analysis from 20 to 60 kg
- LQR-based active suspension control
- Actuator force-saturation analysis
- ISO 8608 Class C random-road simulations
- Vehicle-speed analysis
- Thirty-road statistical evaluation
- Dynamic load coefficient calculation
- Actuator mechanical-power evaluation
- Parameter sensitivity analysis

## MATLAB Files

### `ev_iwm_study.m`

Main MATLAB script used to reproduce the passive and active suspension study, including:

- Centralized-drive, rigid wheel-mounted, and isolated-motor models
- Passive isolation design search
- Fixed and independently retuned motor-mass studies
- Comfort and Balanced LQR controllers
- Bump and Class C random-road responses
- Actuator force-limit study
- Random-road ensemble analysis
- Vehicle-speed analysis
- Parameter sensitivity analysis
- Main manuscript figures and tables

### `ev_force_limited_validation.m`

Standalone MATLAB script for the additional force-limited analyses introduced during peer review, including:

- Dynamic load coefficient (DLC)
- Actuator saturation occurrence
- RMS mechanical power
- Peak absolute mechanical power
- Thirty independent ISO 8608 Class C road realizations
- Force-limited vehicle-speed sweep from 10 to 30 m/s
- Revised force-limited Figure 10

The selected actuator limits are:

- Comfort LQR: ±1500 N
- Balanced LQR: ±1000 N

This script is standalone and does not require the MAT file produced by `ev_iwm_study.m`.

## Requirements

The results were generated and verified using MATLAB R2025b with:

- Control System Toolbox
- Statistics and Machine Learning Toolbox

## Usage

Clone or download this repository and set the repository folder as the current MATLAB folder.

### Main study

Run:

```matlab
ev_iwm_study
```

The script creates:

```text
ev_iwm_data/
ev_iwm_figures/
```

These directories contain the generated numerical results, figures, CSV tables, and MATLAB data file.

### Force-limited validation

Run:

```matlab
ev_force_limited_validation
```

The script creates:

```text
reviewer1_data/
reviewer1_figures/
```

These directories contain the DLC results, saturation and mechanical-power statistics, force-limited speed-sweep data, and revised Figure 10.

## Reproducibility

Random-number seeds used for the passive design search, road-profile generation, road ensemble, and parameter studies are explicitly specified in the MATLAB source code.

The force-limited validation uses:

- Thirty independent Class C road realizations
- A common road phase realization for the vehicle-speed sweep
- A Comfort LQR force limit of ±1500 N
- A Balanced LQR force limit of ±1000 N

Running the scripts with the stated MATLAB version and required toolboxes should reproduce the reported numerical results.

## Citation

Citation information for this repository is provided in `CITATION.cff`.

If this code is used in academic work, please cite the associated manuscript and the archived software release.

## Author

Gridsada Phanomchoeng  
Department of Mechanical Engineering  
Faculty of Engineering  
Chulalongkorn University  
Bangkok, Thailand

## License

This software is released under the MIT License.
