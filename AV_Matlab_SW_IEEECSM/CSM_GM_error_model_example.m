% This file is the companion of the article:
% IMU Error Modeling Tutorial: INS state estimation and real-time sensor calibration 
% authored by 
% Jay A. Farrell, Felipe O. Silva, Farzana Rahman, and J. Wendel.
% This tutorial is accepted for publication in 
% IEEE Control Systems Magazine (referred to below as IEEE CSM)
% The articles main point-of-contact is J. Farrell (farrell@ece.ucr.edu)
%
% This software is distributed for academic purposes under the MIT License
% 
% Copyright (c) 2020 JAY A FARRELL
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is an example of the computations of the Bias instabilty model
% parameters in the paper above. 

clear all
clc
format compact 
format long 

% Method 1
sigma_zB = 7.4e-4;                          % from IEEE CSM Fig. 2, minimum value (almost)
coef     = sqrt((2*log(2))/pi)              % should be 0.664
B_1      = sigma_zB / coef                  % IEEE CSM eqn. (32) solved for B

% Method 2
sigma_zG = 7.4e-4;                          % from IEEE CSM Fig. 2 (flat region)
tau = 60;                                   % from IEEE CSM Fig. 2 (flat region)
T_B = tau/1.89                              % In the same bullet as IEEE CSM eqn. (37)
S_B = (sigma_zG/0.4365)^2/T_B               % IEEE CSM eqn. (37) 
B_2 = sqrt(S_B*pi*(0.4365)^2*T_B/(2*log(2))) % solving IEEE CSM eqn. (39) for B

