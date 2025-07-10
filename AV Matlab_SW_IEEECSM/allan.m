function [avar]=allan(data,tau)
% Compute various Allan deviations for a constant-rate time series
% [avar]=allan(data,tau) 
%
% INPUTS:
% data should be a struct containing the following fields: 
%  data.freq    the time series measurements in arbitrary units
%  data.rate    constant rate of time series in (Hz)
%               (Differently from previous versions of allan.m,
%               it is not possible to compute variances for time-
%               stamp data anymore.)
% tau is an array of the tau values for computing Allan deviations
%
% OUTPUTS: 
% avar is a struct containing the following fields (for values of tau):
%  avar.sig     = standard deviation
%  avar.sig2    = Allan deviation
%  avar.sig2err = standard error of Allan deviation
%  avar.osig    = Allan deviation with overlapping estimate
%  avar.osigerr = standard error of overlapping Allan deviation
%  avar.msig    = modified Allan deviation 
%  avar.msigerr = standard error of modified Allan deviation
%  avar.tsig    = timed Allan deviation
%  avar.tsigerr = standard error of timed Allan deviation
%  avar.tau1    = measurement interval in (s)
%  avar.tauerr  = errors in tau that might occur because of initial
%  rounding
%
% NOTES:
% Calculations of modified and timed Allan deviations for very long time
% series become very slow. It is advisable to uncomment .msig* and .tsig*
% only after calculations of .sig*, .sig2* and .osig* have been proven
% sufficiently fast.
%
% No pre-processing of the data is performed.
% For constant-rate time series, the deviations are only calculated for tau
% values greater than the minimum time between samples and less than half
% the total time.
%
% versionstr = 'allan v3.0';
% FCz OCT2009
% v3.0  faster and very plain code, no plotting; various Allan deviations
%       can be calculated; script and sample data are availabie on
%       www.nbi.dk/~czerwin/files/allan.zip
%		(Normal, overlapping and modified Allan deviations are calculated 
%       in one function, in strong contrast to MAHs approach of splitting 
%       up among various functions. This might be beneficial for individual 
%       cases though.)
% 
% MAH 2009
% v2.0 and others
%
% FCz OCT2008
% v1.71 'lookfor' gives now useful comments; script and sample data are
%       availabie on www.nbi.dk/~czerwin/files/allan.zip
% v1.7  Improve program performance by mainly predefining matrices outside
%       of loops (avoiding memory allocation within loops); no changes to
%       manual
%
% early program core by Alaa MAKDISSI 2003
% (documentation might be found http://www.alamath.com/)
% revision and modification by Fabian CZERWINSKI 2009
%
% For more information, see:
% [1] Fabian Czerwinski, Andrew C. Richardson, and Lene B. Oddershede,
% "Quantifying Noise in Optical Tweezers by Allan Variance,"
% Opt. Express 17, 13255-13269 (2009)
% http://dx.doi.org/10.1364/OE.17.013255

n=length(data.freq);
jj=length(tau);
m=floor(tau*data.rate);

avar.sig     = zeros(1, jj);
avar.sigerr  = zeros(1, jj);
avar.sig2    = zeros(1, jj);
avar.sig2err = zeros(1, jj);
avar.osig    = zeros(1, jj);
avar.osigerr = zeros(1, jj);

tic;

for j=1:jj
        
    D=zeros(1,n-m(j)+1);
    D(1)=sum(data.freq(1:m(j)))/m(j);
    for i=2:n-m(j)+1
        D(i)=D(i-1)+(data.freq(i+m(j)-1)-data.freq(i-1))/m(j);
    end
    
    %standard deviation
    avar.sig(j)=std(D(1:m(j):n-m(j)+1));
    avar.sigerr(j)=avar.sig(j)/sqrt(n/m(j));
    
    %normal Allan deviation 
    avar.sig2(j)=sqrt(0.5*mean((diff(D(1:m(j):n-m(j)+1)).^2)));
    avar.sig2err(j)=avar.sig2(j)/sqrt(n/m(j));
    
    %overlapping Allan deviation
    z1=D(m(j)+1:n+1-m(j));
    z2=D(1:n+1-2*m(j));
    u=sum((z1-z2).^2);
    avar.osig(j)=sqrt(u/(n+1-2*m(j))/2);
    avar.osigerr(j)=avar.osig(j)/sqrt(n-m(j));
    
    % toc
    
end

avar.tau1=m/data.rate;
avar.tauerr=tau-avar.tau1;

toc;
end