%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD MIP DSA Functor
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Generates the LHS MIP diffusion matrix, A, and the RHS
%                   scattering matrix, S.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,S] = LD_MIP_DSA_func(lambda, data)
% Get Geometric and XS Info
dx = data.dx;
M = data.mats.M;
K = data.mats.K;
eL = data.mats.eL; eR = data.mats.eR;
bL = data.mats.bL; bR = data.mats.bR;
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
kpen = max(c*pp*(pp+1)*(D/dx + D/dx)/2,.25);
nL = -1; nR = 1;
S = (sigs*M) * PM;
A = (siga*M + D*K) * PM;
% ( [[u]] , [[b]] )
A = A + kpen*[1,-1;-e,e];
% ( {{Du}} , [[b]] )
A = A + nL*D/2*[-(e0+em);(e0+e)]/dx*(-bL);  % left boundary
A = A + nR*D/2*[-(e+e0);(ep+e)]/dx*(-bR);   % right boundary
% ( [[u]] , {{Db}} )
A = A + nL*D/2*eR*[e0,-e0];        % left boundary
A = A + nR*D/2*eL*[e,-e];          % right boundary