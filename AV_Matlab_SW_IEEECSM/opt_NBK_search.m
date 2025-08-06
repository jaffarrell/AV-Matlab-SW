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

% The following functions implement the optimization-based approach
% described in subsection entitled "Optimization-based Parameter Selection
% for a Given Model"

function []= opt_NBK_search()
format long
format compact
clc

COMPUTE_ASD = 0;      % only compute the ASD when tau changes; otherwise load it.
if COMPUTE_ASD == 1
    % NAV1000 IMU accelerometer z axis
    load parsed_isolated_marble_data_az
    L = length(acc_z);        % Length of the data set.
    Fs = 100;                 % IMU sampling rate = 100Hz
    T = 1/Fs;                 % time difference
    
    %% setup to compute AV and ASD
    data.freq = acc_z;        % computing allan variance plot
    data.rate = Fs;
    N_tau = 50;             % number of tau values at which to optimize
    tau_min = 0.04;         % min value of tau range
    tau_max = 500;          % max value of tau range
    tau = logspace(log10(tau_min/T),log10(tau_max/T),N_tau)*T;
    % compute ASD
    sprintf('Computing ASD')
    [avar] = allan(data,tau);
    asd(1:N_tau,1) = avar.sig2;
    [m,n] = size(asd);
    if n>m
        asd = asd'; % ensure it is a column
    end
    av(1:N_tau,1) = (asd.^2);
    % Compute weighting matrix for cost function
    er  = avar.sig2err;                     % ASD standard deviation
    c = 4 * avar.sig2err .* avar.sig2err;   % Discussion near IEEE CSM eqn. (15)
    W = inv(diag(c));                       % convert to diagonal weight matrix
    sprintf('Finished computing ASD')
    save 'ASD_Marble_slab.mat' tau asd av W er c
elseif COMPUTE_ASD == 2  % Use data that can be perfectly fit (for debugging)
    sprintf('Computing ASD: known case')
    N = 3.3e-3;                                % random walk parameter [m/sec^(3/2)]
    K = 1.2e-4;                                % rate random walk parameter [m/sec^(5/2)]
    B = 6e-4;                                  % Bias instability parameter [m/sec^2]
    TB = 50;
    
    SN = N^2;                                  % IEEE CSM eq. (22)
    SK = K^2;                                  % IEEE CSM eq. (25)
    SB = (2*B^2*log(2))/(pi*0.4365^2*TB);      % IEEE CSM eq. (39)
    
    L = 360000;      % artificial data length. Only affects error bound.
    
    N_tau = 50;
    tau_min = 0.04;
    tau_max = 500;
    Fs = 100;                                  % IMU sampling rate = 100Hz
    T = 1/Fs;
    tau_vec = logspace(log10(tau_min/T),log10(tau_max/T),N_tau)*T;
    % Define vector of cluster times
    
    % Generate Allan Variance plot data
    av = zeros(N_tau,1);    % preallocate
    for k=1:N_tau
        tau = tau_vec(k);
        av(k,1) = SB*TB^2/tau * (1 - TB/(2*tau)*(3 - 4*exp(-tau/TB) + exp(-2*tau/TB))) + SN/tau + SK*tau/3;
        asd(k,1) = sqrt(av(k,1));
    end
    % following are defined for compatibility with following code
    tau = tau_vec;        %
    v(1:N_tau,1) = 1./sqrt( 2*(L*T./tau -1) ); % percentage uncertainty dfnd in IEEE CSM eqn. (15)
    er= v.*asd;           % absolute uncertainty (treat this as standard dev)
    c = (er.*er);         % absolute variance
    W = inv(diag(c));     % convert to diagonal weight matrix
else
    sprintf('Loading ASD from file')
    load ASD_Marble_slab.mat
end %if

%% Manually extracted parameters from Table 2 in IEEE CSM 
T_B = 20;                % correlation time [sec]
% mu_B = 1/T_B;          % time constant parameter [Hz]
N = 3.3e-3;              % random walk parameter [m/sec^(3/2)]
K = 1.4e-4;              % rate random walk parameter [m/sec^(5/2)]
B = 4e-4;                % Bias instability parameter [m/sec^2]

S_N = N^2;                                  % IEEE CSM eq. (22)
S_K = K^2;                                  % IEEE CSM eq. (25)
S_B = (2*(B^2)*log(2))/(pi*0.4365^2*T_B);   % IEEE CSM eq. (35)
theta_o  = [S_B,T_B,S_N,S_K]';              % IEEE CSM eq. (63) with different order
phi      = cmpt_basis(tau,theta_o(2));
[av_o]   = cmpt_AV(theta_o,phi);
[Cost_o] = cmpt_cost(av,av_o,W);

%% Optimal parameter search
%    Golden section search over theta(2). See Wikopedia for description.
%    To find the minimum Cost on [a,b]
%

% range of T_B in [a,b] seconds
a=10;
[theta_a,Cost_a]=cmpt_theta(tau,av,a,W);
b=200;
[theta_b,Cost_b]=cmpt_theta(tau,av,b,W);
mu =  (sqrt(5) + 1)/2;      % Golden ratio
tol=  0.01;
cnt = 0;
while abs(b - a) > tol
    c = b - (b - a) / mu;
    d = a + (b - a) / mu;
    [theta_c,Cost_c]=cmpt_theta(tau,av,c,W);
    [theta_d,Cost_d]=cmpt_theta(tau,av,d,W);
    [a,c,d,b
    Cost_a, Cost_c, Cost_d, Cost_b] % display search process
    % implement golden section search
    if Cost_c < Cost_d
        b = d;
        Cost_b = Cost_d;
        plot_AV_mdl_components(tau,av,theta_c,er)
    else
        a = c;
        Cost_a = Cost_c
        plot_AV_mdl_components(tau,av,theta_d,er)
    end
    pause_time = 1.0;   % seconds
    pause(pause_time)
end
[a,b
 Cost_a, Cost_b] % display search process
theta_v = theta_c;                  % computed optimal value
phi      = cmpt_basis(tau,theta_v(2));
[av_v]=cmpt_AV(theta_v,phi);        % optimized fit to AV


AV_SBSNSK = [theta_v([1,3,4]),theta_o([1,3,4])]              % display result
ASD_NK = sqrt([theta_v([3,4]),theta_o([3,4])])              % display result
TB = [theta_v(2),theta_o(2)]

plot_AV_error(tau,av,av_v,av_o,W)

asd_o = sqrt(av_o);
asd_v = sqrt(av_v);

figure(1);clf
asd_min = 1e-4;
asd_max = 0.1;
loglog(tau,asd,'b-','LineWidth',2)
hold on
grid on
loglog(tau,asd_o,'g--','LineWidth',2)
loglog(tau,asd_v,'r--','LineWidth',2)
xlabel('cluster time (sec)','FontSize',14)
ylabel('ASD (m/s^2)','FontSize',14)

lower= max([asd-er,0.1*asd_min*ones(size(asd))]')';
upper= min([asd+er, 10*asd_max*ones(size(asd))]')';
loglog(tau,upper,'b.',tau,lower,'b.','LineWidth',2)

ylim([asd_min asd_max])
legend('Data','Paper Model','Optimal Model (min AVAR sq err)')


function []=plot_AV_error(tau,av_t,av_c,av_o,W)
% plot the AV and errors
% av_t --- computed from data
% av_c --- optimal fit by computation
% av_o --- from paper
err =  (av_t - av_c);   % curve fit error raw
figure(3)
clf
subplot(311)
loglog(tau,av_t,'b-','LineWidth',2)
grid on
hold on
loglog(tau,av_c,'r--',tau,av_o,'g--','LineWidth',2)
xlabel('cluster time (sec)','FontSize',14)
ylabel('AVAR','FontSize',14)
%ylim([1e-8 1e-2])
legend('Data','Optimized Mdl','Paper Mdl')  % Optimize based n (min AVAR sq err)

subplot(312)
loglog(tau,abs(err),'r.',tau,abs(av_t-av_o),'g.','LineWidth',2)
grid on
xlabel('cluster time (sec)','FontSize',14)
ylabel('AVAR raw Error','FontSize',14)
legend('curve fit','paper')


subplot(313)
loglog(tau,W*abs(err),'r.',tau,W*abs(av_t-av_o),'g.','LineWidth',2)
grid on
xlabel('cluster time (sec)','FontSize',14)
ylabel('AVAR weighted Error','FontSize',14)
legend('curve fit','paper')

function []=plot_AV_mdl_components(tau,av,theta,err)
TB = theta(2);
[phi]=cmpt_basis(tau,TB);
[av_T]=cmpt_AV(theta,phi);
av_B = phi(:,1)*theta(1);
av_N = phi(:,3)*theta(3);
av_K = phi(:,4)*theta(4);
figure(4); clf
loglog(tau,av,'.')
xlabel('Delay, \tau, s')
ylabel('AV')
grid
ylim([1e-8 1e-4])
hold on
loglog(tau,[av_N,av_B,av_K,av_T])
hold off
legend('AV actual','AV_N','AV_B','AB_K','AV-fit')
title('AV optimization-based fitting process')

figure(5)
loglog(tau,sqrt(av),'.')
xlabel('Delay, \tau, s')
ylabel('ASD')
grid
ylim([1e-4 1e-2])
hold on
loglog(tau,sqrt([av_N,av_B,av_K,av_T]))
hold off
legend('ASD actual','ASD_N','ASD_B','ASD_K','ASD-fit')
title('ASD optimization-based fitting process')

figure(15)
errorbar(tau,sqrt(av),err);
ylim([1e-4 1e-2])
set(gca,'XScale','log','YScale','log');
xlabel('Delay, \tau, s')
ylabel('ASD')
grid
hold on
loglog(tau,sqrt([av_N,av_B,av_K,av_T]))
hold off
legend('ASD actual','ASD_N','ASD_B','ASD_K','ASD-fit')
title('ASD optimization-based fitting process')

function [Cost]=cmpt_cost(av_t,av_c,W)
% compute the cost functions that we are trying to minimize
Cost =  (av_t - av_c)'*W*W*(av_t - av_c);

function [phi]=cmpt_basis(tau,theta_2)
% follows notation preceeding IEEE CSM eqn (63) in a different order:
%   theta(1) is IEEE CSM theta_2
%   theta(2) is IEEE CSM theta_4 (nonlinear in IEEE CSM eqn. (63))
%   theta(3) is IEEE CSM theta_1 
%   theta(4) is IEEE CSM theta_3 
% phi(:,2) is not needed here because theta(2) is not linear

a(:,1) = (theta_2./tau)';       % worker variable
b      = ones(size(a));         % worker variable
c      = b./a;                  % worker variable
phi(:,1) = theta_2*a.*(b-(a/2).*(3*b-4*exp(-c)+exp(-2*c)));
phi(:,3) = 1./tau;      % coefficeint of theta(3)
phi(:,4) = tau./3;      % coefficient of theta(4)

function [theta,Cost]=cmpt_theta(tau,av,theta_2,W)
% Given a value for theta(2), optimize the other parameters that appear linearly

[phi]      = cmpt_basis(tau,theta_2);
H          = phi(:,[1,3,4]);                % discard column 2
lambda     = inv(H'*W*H)*H'*W*av;           % LS over theta 1, theta 3, theta 4
theta      = [lambda(1);theta_2;lambda(2);lambda(3)];
[avc]      = cmpt_AV(theta,phi);
[Cost]     = cmpt_cost(av,avc,W);


function [av]=cmpt_AV(theta,phi)
% follows notation preceeding IEEE CSM eqn (63) in a different order:
%   theta(1) is IEEE CSM theta_2
%   theta(2) is IEEE CSM theta_4 (nonlinear in IEEE CSM eqn. (63))
%   theta(3) is IEEE CSM theta_1 
%   theta(4) is IEEE CSM theta_3 
% phi(:,2) is zero
% theta(2) is accounted for in phi(1)
lambda    = [theta(1);theta(3);theta(4)];
av(:,1)  = phi(:,[1,3,4])*lambda;