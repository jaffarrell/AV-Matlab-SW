This directory contains the following m-files:

-- allan.m: Implements the Allan variance computation as a function.  The authors of this function are documented in the header to the file. 

-- ASD_plot_generation.m: Given parameters that define a state-space model for the sensor stochastic errors, this function: 
(1) computes the state-space model in continuous-time and converts to discrete-time, 
(2) simulates that model to produce a stochastic error sequence, 
(3) computes and plots the ASD.

-- CSM_GM_error_model_example.m: computes B from sigma_B and sigma_G

-- opt_NBK_search.m: Implements an optimization based approach to selecting the state-space model parameters for the problem defined in Section VIII.B of the paper. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


The following data files are important:

The only author supplied data file is:
-- parsed_isolated_marble_data_az.mat: This is data from a motion isolated IMU. It is used by "opt_NBK_search.m". 

The following data files are computed by the software and saved to speed-up repeated computations:
-- ASD_Marble_slab.mat: This data file is output by "opt_NBK_search.m" to save repeated computation of the ASD, which is slow. The m-file contains a flag "COMPUTE_ASD". This flag should be set to 1 the first time the m-file is run to create this file. After that, it can be reset to 0 to save computations.  
 
-- data_Z.mat: This file is both created and read by ASD_plot_generation.m.  This m-file contains a flag "rerun". Set the flag to 1 the first time that you run the code, and any time thereafter when you want new data, to simulate the model and create the file.  Once you have the file, you can reset the flag to 0 to save computation time. Note that if you change the state-space error model, you need to set the flag to generate new data for the new model. 

