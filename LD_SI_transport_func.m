%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD MIP DSA Functor
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Generates the flux moment LHS system matrix, T.
%
%                   T = D*L^(-1)*M*S
%
%                   D = Discrete-to-Moment Operator
%                   L = Transport Operator (streaming + interaction)
%                   M = Moment-to-Discrete Operator
%                   S = Scattering Operator
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = LD_SI_transport_func(lambda, data)
% Get Geometric and XS Info
dx = data.dx;
M = data.mats.M;
G = data.mats.G;
sigt = data.XS.sigt;
sigs = data.XS.sigs;
% Build Phase Info
e = exp(1i*lambda*dx);
PM = [1,0;0,e];
% Build All Transport/Diffusion Matrices
% --------------------------------------
S = (sigs*M/data.SN.norm) * PM;
T = zeros(data.ndofs);
% Loop through Quadrature and Build Transport Matrix
for m=1:data.SN.num_dirs
    qw = data.SN.w(m);
    qx = data.SN.mu(m);
    % Calculate Transport Volume Contribution
    L = (sigt*M - qx*G) * PM;
    % Add Transport Edge Contributions
    if qx > 0
        L = L + qx*[0,-1;0,e];
    else
        L = L + qx*[-1,0;e,0];
    end
    % Add Angular Transport Contribution
    T = T + qw*(L\S);
end