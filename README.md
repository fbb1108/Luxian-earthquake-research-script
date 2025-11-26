This study employed Wang Rongjiang’s POEL software (https://github.com/RongjiangWang/POEL_2024) to perform the poroelastic simulations. The resulting snapshot data were used to compute stress and pore pressure variations in cylindrical coordinates.

The provided MATLAB scripts (c1–c5) sequentially convert the cylindrical-coordinate stress fields into Cartesian coordinates, enabling the calculation of Coulomb stress changes on the target fault plane.

The script combine_coulomb_polar2NED_opt.m is designed to superimpose the effects of multiple injection platforms and to determine the optimal fault-slip mechanism under combined stress perturbations.

All scripts were developed for post-processing POEL simulation outputs and generating the stress and fault response results presented in the manuscript.
