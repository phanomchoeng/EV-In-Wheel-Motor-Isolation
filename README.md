# EV In-Wheel Motor Isolation

MATLAB source code accompanying the manuscript:

**"Mitigating Unsprung Mass in an Electric Vehicle through Passive and Active Isolation of a Wheel-Integrated Motor"**

## Description

This repository contains the MATLAB source code used to reproduce the numerical analyses presented in the manuscript.

The code includes:

- Mass-matched centralized-drive, rigid in-wheel-motor, and isolated-motor EV models
- Passive motor-isolation parameter search
- Motor-mass analysis from 20 to 60 kg
- LQR-based active suspension control
- Actuator force-saturation analysis
- ISO 8608 Class C random-road simulations
- Vehicle-speed analysis
- Random-road ensemble analysis
- Parameter sensitivity analysis

## Usage

Run the following MATLAB script:

    ev_iwm_study.m

The script automatically creates the required output directories and generates the numerical results, figures, CSV tables, and MATLAB data file.

## Requirements

MATLAB with:
- Control System Toolbox
- Statistics and Machine Learning Toolbox

## Reproducibility

Random seeds used in the numerical analyses are explicitly specified in the source code to support reproducibility of the reported results.

## Author

Gridsada Phanomchoeng  
Department of Mechanical Engineering  
Faculty of Engineering  
Chulalongkorn University  
Bangkok, Thailand

## License

This software is released under the MIT License.
