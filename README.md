# NMMD / MDSPACE --  *Extension of GENESIS 2.1*

## Installation: 
*See https://www.r-ccs.riken.jp/labs/cbrt/installation/ for installation requirements*

To build GENESIS :
```
cd /path/to/genesis
autoreconf -fi
./configure
make install
```

## Usage NMMD:

The input file is similar as in GENESIS (see how to define the input file in https://www.r-ccs.riken.jp/labs/cbrt/usage/). To use NMMD, you must chose the `integrator=NMMD` in the `DYNAMICS` section and add a section `NMMD` to the INP file as follow :

```
...

[DYNAMICS]
integrator = NMMD

...

[NMMD]
nm_number = 10
nm_mass = 10.0
nm_file = /path/to/nm/file
nm_dt = 0.001
```
- `nm_number` : The number of normal modes to use after skipping the 6 first modes (example : 10 will use modes from 7-16), default 10
- `nm_mass` : Mass value of Normal modes for NMMD, default : 10.0
- `nm_file` : file containing the normal mode vectors
- `nm_dt` : normal mode integration time step, default :0.001
 

## Limitations:
- NMMD is available only for ATDYN
- NMMD is available only for LANGEVIN temperature control in the NVT enemble
- SHAKE/RATTLE algorithm have to be turn off

## Usage MDSPACE:
The MDSPACE extension allows to fit images instead of volumes

usage : 

```
[EXPERIMENTS]
emfit = YES                            # YES/NO
emfit_type = IMAGE                     # VOLUME/IMAGE 
emfit_target = path/to/image/file.spi  # 2D SPIDER file
emfit_sigma = 2.0                      # sigma of 2D gaussian
emfit_tolerance = 0.01                 # Gaussian truncation threshold
emfit_period = 1                       # not used
emfit_roll_angle = 0.0                 # Euler roll angle in degrees
emfit_tilt_angle = 0.0                 # Euler tilt angle in degrees
emfit_yaw_angle = 0.0                  # Euler yaw angle in degrees
emfit_shift_x = 0.0                    # Shift in x direction in pixels
emfit_shift_y = 0.0                    # Shift in y direction in pixels
emfit_pixel_size = 1.0                 # Size of a pixel in Angstrom
```
## Limitations:
- NMMD is available only for ATDYN

# genesis-2.0
Source code for genesis 2.0 beta for Fugaku
For GENESIS compile on Fugaku and execution, please see READMD.howto
