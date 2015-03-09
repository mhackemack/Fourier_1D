%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get 1D Local Matrices
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Returns the FE matrices by a given cell width (dx).
%
%                   M = Cell Mass Matrix
%                   K = Cell Stiffness Matrix
%                   G = Cell Gradient Matrix
%                   eL = Left-Face Gradient matrix
%                   eR = Right-Face Gradient matrix
%                   bL = Left-Face Mass matrix
%                   bR = Right-Face Mass matrix
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mats = get_1D_mats(dx)
mats.M = [2,1;1,2]*dx/6;
mats.K = [1,-1;-1,1]/dx;
mats.G = [-1,-1;1,1]/2;
mats.eL = [-1;1]/dx;
mats.eR = [1;-1]/dx;
mats.bL = [1,0];
mats.bR = [0,1];