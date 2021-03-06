%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD IP DSA Functor
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Generates the LHS IP diffusion matrix, A, and the RHS
%                   scattering matrix, S.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,S] = LD_IP_DSA_func(lambda, data)
% Get Geometric and XS Info
dx = data.dx;
M = data.mats.M;
K = data.mats.K;
bL = data.mats.bL;
bR = data.mats.bR;
siga = data.XS.siga;
sigs = data.XS.sigs;
D = data.XS.D;
% Build Phase Info
e = exp(1i*lambda*dx);
e0 = 1; em = e0/e; ep = e*e;
PM = [1,0;0,e];
% Calculate Diffusion Matrices
% ----------------------------
pp = 1; c = data.IPConstant;
kpen = c*pp*(pp+1)*(D/dx + D/dx)/2;
nL = -1; nR = 1;
S = (sigs*M) * PM;
A = (siga*M + D*K) * PM;
% ( [[u]] , [[b]] )
A = A + kpen*[1,-1;-e,e];
% ( {{Du}} , [[b]] )
A = A + nR*D/2*[-(e+e0);(ep+e)]/dx*(-bR);   % right boundary
A = A + nL*D/2*[-(e0+em);(e0+e)]/dx*(-bL);  % left boundary
% ( [[u]] , {{Db}} )
A = A + nL*D/2*eR*[e0,-e0];        % left boundary
A = A + nR*D/2*eL*[e,-e];          % right boundary