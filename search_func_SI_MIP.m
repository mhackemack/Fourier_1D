%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD SI+MIP minsearch functor
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Functor that returns the negative of the largest eigenvalue 
%                   for the SI+MIP system. This acts like a minimizing function.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = search_func_SI_MIP(lambda, data)
% Retrieve Matrices
T = LD_SI_transport_func(lambda, data);
[A, B] = LD_MIP_DSA_func(lambda, data);
% Calculate eigenspectrum and output min results
I = eye(data.ndofs);
P = T + (A\B)*(T-I);
out = - max(abs(eig(P)));