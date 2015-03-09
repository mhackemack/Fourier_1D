%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Get Logarithmic Spacing
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    Get a uniform spacing in log-space.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        If a 2nd output is specified, the routine provides the
%                   tick-mark locations for the major divisions. This was
%                   necessary for pretty-plot pictures. If this second output is
%                   desired, make sure your min/max x-values lie on integer
%                   values in log-space (e.g. 1e-4, 1e0, 1e2, etc.). I'm not
%                   sure how it will make the plots look otherwise.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_logarithmic_x( xmin, xmax, num )
lmin = log10(xmin);
lmax = log10(xmax);
lspace = linspace(lmin, lmax, num);
out = 10.^(lspace);
% Set Outputs
varargout{1} = out;
if nargout == 2
    xticks_log = round(lmin):round(lmax);
    varargout{2} = 10.^(xticks_log);
end