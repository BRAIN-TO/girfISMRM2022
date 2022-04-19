# Demonstration Code for Gradient Impulse Response Function (GIRF) Calculation and Analysis

This is the demonstration code of the abstract *MR Gradient System Long-Term Stability Investigation and Protocol Optimization for Quality Control using Gradient Impulse Response Function (GIRF)*`* for the Joint Annual Meeting ISMRM-ESMRMB & ISMRT 31st Annual Meeting, which will be held in London UK, between May 7 - 12, 2022. The program number and the link will be added once it is announced.

> Copyright (C) 2021-2022  
> Zhe "Tim" Wu  
> <zhe.wu@rmp.uhn.ca>  
>  
> Brain Research in Advanced Imaging and Neuromodeling - Toronto (BRAIN-TO) Lab  
> Techna Institute
> University Health Network

## Quick Start

1. Download the demonstration data from [here](https://www.doi.org/10.5281/zenodo.6376737) and extract them under the same folder. See details of the description of each file in the section [Dataset for Demonstration](#dataset-for-demonstration);

2. In the `main.m` file, change the variables under the section "Set User Parameters". Set the variable `dataPath` as the path that stores all the downloaded data; set the variable `gradientAxis` to the gradient axis that needs to be analysed (choose from 'X', 'Y', and 'Z', **case insensitive**); set the variable `measNum` to the measurement dataset to be analysed (only applicable to `script_OriginGIRFCalculation.m`, see the following section [Script List](#script-list) and the instructions inside this script for details.)

3. Run the `main.m` file and all the demo scripts should be executed. You may change the variable `gradientAxis` to check the other gradient axis. 

## Dataset for Demonstration

The demonstration data used by the code in this repository can be found [here](https://www.doi.org/10.5281/zenodo.6376737). **There are three .zip files and they should be copied to and be extracted under the same folder.**

The demonstration data contains three parts:

1. `Meas1.zip` contains raw T2* decay signal acquired only with the positive triangle blips using phantom-based method (See Figure 1 in the abstract).

2. `Meas2.zip` contains raw T2* decay signal acquired with both positive and negative triangle blips using phantom-based method （See Figure 1 in the abstract）. This data was acquired 8 months after the acquisition of `Meas1.zip`.

3. `CalculatedGIRF.zip` provides the author's pre-calculated GIRFs using the data without coil averaging. This data is used for all the postprocessing (e.g. SNR and stability analysis, comparison between original and optimized GIRF methods, etc.) in the published abstract with the source code provided in the same Github repository.

Please note that the raw T2* decay signal in `Meas1.zip` and `Meas2.zip` is only for the purpose of demonstrating GIRF calculations (using blips with single polarity and dual polarities) - they are NOT used in the analysis in the abstract. The reason is that we averaged across the coil dimension of raw to save data volume for demonstration purposes, which leads a reduced SNR of the calculated output gradients and GIRFs. For the analysis of GIRF (SNR, stability, etc.) we use the pre-calculated GIRFs in `CalculatedGIRF.zip`.

## Script List

1. `main.m` executes all the following scripts. The users need to set the path of the downloaded data and the axis they would like to calculate or analysis. Alternatively, each following scripts can be run individually.

2. `script_OriginGIRFCalculation.m` demonstrates the calculation of GIRF using the original method (positive blips only).

3. `script_OptimizedGIRFCalculation.m` demonstrates the calculation of GIRF using the optimized method (both positive and negative blips).

4. `script_GIRFStability.m` analyses stability of GIRFs in 8-month gap (Figure 3 C-E).

5. `script_GIRFDiffBetweenMethods.m` compares the GIRF results from two calculation methods (Figure 4 first row).

6. `script_SNRAnalysis.m` analyses the SNR level of the GIRF results from two calculation methods (Figure 4 first row).

## Execution Efficiency and Required Resources

The code in this repository were tested on a laptop equipped with Intel Core i7-8850H @ 2.60 GHz.

Execution time (all scripts through `main.m`): < 10 seconds.

Peak RAM capacities needed: < 3GB.

## References

**Gradient system characterization**
Vannesjo, S.J., Haeberlin, M., Kasper, L., Pavan, M., Wilm, B.J., Barmet, C., Pruessmann, K.P., 2013. Gradient system characterization by impulse response measurements with a dynamic field camera. Magn Reson Med 69, 583–593. https://doi.org/10.1002/mrm.24263

**Image reconstruction based on the GIRF characterization**
Vannesjo, S.J., Graedel, N.N., Kasper, L., Gross, S., Busch, J., Haeberlin, M., Barmet, C., Pruessmann, K.P., 2016. Image reconstruction using a gradient impulse response model for trajectory prediction. Magn Reson Med 76, 45–58. https://doi.org/10.1002/mrm.25841

**GIRF-based spiral fMRI**
Graedel, N.N., Kasper, L., Engel, M., Nussbaum, J., Wilm, B.J., Pruessmann, K.P., Vannesjo, S.J., 2019. Feasibility of spiral fMRI based on an LTI gradient model. bioRxiv 805580. https://doi.org/10.1101/805580

**GIRF-based pre-emphasis of gradient or shim channels**
Vannesjo, S.J., Duerst, Y., Vionnet, L., Dietrich, B.E., Pavan, M., Gross, S., Barmet, C., Pruessmann, K.P., 2017. Gradient and shim pre-emphasis by inversion of a linear time-invariant system model. Magn Reson Med 78, 1607–1622. https://doi.org/10.1002/mrm.26531
