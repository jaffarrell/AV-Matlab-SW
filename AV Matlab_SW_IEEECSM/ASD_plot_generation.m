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

%% This script completes the following tasks:
% (1) Given a continuous time error model, it computes a 
% discrete-time equivalent state space-model.
% (2) Given that discrete-time model, it produces a
% stochastic error sequence.
% (3) Given a sequence of stochastic errors, computes and
% plots the Allan Variance data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%%-------------------------%
% Task 1 - Define various parameters, set-up the continuous-time model, and
%          convert it to a discrete-time model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sensor sampling frequency
Fs=100;                 % sample frequency, Hz (instrument dependent)
dT=1/Fs;                % sample time, sec

% Data length
L_n=10000000;   % number of samples (experiment dependent)

% IMU time stamp
imu_t=0:dT:(L_n-1)*dT;

% Random walk parameter, N (m/sec.sqrt(sec))
% See IEEE CSM section containing eqns. (21-23) 
% See IEEE CSM Table 2 
N=3.3e-3;               % (instrument dependent)
% PSD parameter for random walk process noise
% See IEEE CSM eqn. (22)
Sn=N^2;

% Rate random walk parameter, K (m/sec.sec^(3/2))
% See IEEE CSM section containing eqns. (24-27) 
% See IEEE CSM Table 2 
K=1.4e-4;               % (instrument dependent)
% PSD parameter for rate random walk
% See IEEE CSM eqn. (25)
Sk=K^2;

% Bias instability parameter, B (m/sec.sec)
% See IEEE CSM section containing eqns. (31-32) 
% See IEEE CSM Table 2 
B=0.0005;               % (instrument dependent)

% Correlation time in sec for the Gauss-Markov Error Model of z_G(t)
% See IEEE CSM Table 2 
Tc=32;                  % T_B in eqn. (34) (instrument dependent)
% Gauss-Markov pole location is as s = -mu_B
mu_b=1/Tc;              % IEEE CSM eqn. (34) 

% Continuous-time PSD parameters
% See IEEE CSM eqn. (39)
% PSD parameter for bias instability
Sb=(2*B^2*log(2))/(pi*0.4365^2*Tc);

% Continuous-time state-space model 
% See IEEE CSM eqn. (40)
A_t=[-mu_b 0;
    0 0];
B_t=[1 0;
    0 1];
C_t=[1 1];

% Process noise for noise propagation 
% See IEEE CSM eqn. (41) 
S_w_z=[Sb 0;
        0 Sk];

% Computation of discrete-time equivalent stae-space model
%
% [3] refers to IEEE CSM reference [3]: J. A. Farrell, "Aided Navigation:
% GPS with High Rate Sensors, McGraw-Hill, 2008.
%
% Continuous to discrete transformation using eqn. (4.114) of [3]
A_c(1:2,1:2) = -1*A_t;
A_c(1:2,3:4) = B_t * S_w_z * B_t';
A_c(3:4,3:4) = A_t';

A_c = A_c*dT;
B_c = expm(A_c);

% From eqn. (4.115) of [3]
Phi = B_c(3:4,3:4)';

% From eqn. (4.116) of [3]
Qd  = Phi * B_c(1:2,3:4);

% Discrete time bias instability driving noise
Qbk = Qd(1,1);

% Discsrete time rate random walk process noise
Qkk = Qd(2,2);

% Discsrete time random walk process noise
% See IEEE CSM eqn. (62)
Qnk = Sn*Fs;

%-------------------------%
% Task 2 -Simulate the discrete-time system
% There are two modes to run: 
% rerun = 0 (Load the previously saved generated data)
% rerun = 1 (Generate data using given parameters)
% Set rerun to 1 the first time that you run the code. 
% Set rerun to 0 for subsequent runs, unless you want a new data sequence.
rerun = 1;
if rerun == 1
    display('Simulating discrete-time model')
    sQbk = sqrt(Qbk);
    sQkk = sqrt(Qkk);
    sQnk = sqrt(Qnk);      
    Z    = zeros(L_n,1);    % preallocation
    % Initial error state
    Z_b=[0;0];
    for ij=1:L_n
        % Generate random noise 
        n_b = sQbk * randn(1,1);
        n_k = sQkk * randn(1,1);
        n_n = sQnk * randn(1,1);
        % IEEE CSM eqn. (44) 
        Z_b = Phi*Z_b + [n_b;n_k];
        % IEEE CSM eqn. (45) 
        Z(ij,1)= C_t * Z_b + n_n;
    end
    save data_Z.mat Z
else
    load data_Z.mat
end

%-------------------------%
%--------Task 3-----------%
display('Computing ASD for simulated data')
% Computing Allan variance plot
data.freq=Z;
data.rate=Fs;
% Define the cluster times
tau1=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09];
tau2=10*dT:10*dT:90*dT;
tau3=100*dT:100*dT:900*dT;
tau4=1000*dT:1000*dT:9000*dT;
tau5=10000*dT:10000*dT:90000*dT;
tau6=100000*dT:100000*dT:900000*dT;
tau7=1000000*dT:1000000*dT:4000000*dT;
tau=[tau1 tau2 tau3 tau4 tau5 tau6 tau7]';
% Compute AV
[avar]=allan(data, tau);

% Allan SD is in sig2 field of avar structure
allan_sd=avar.sig2;
figure
loglog(tau,allan_sd,'b.','LineWidth',2)
xlabel('Cluster time (sec)','FontSize',14)
ylabel('Allan SD (m/s^2)','FontSize',14)
ylim([1e-4 1e-1]);
grid on;
title('ASD plot from generated data')





