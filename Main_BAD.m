%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD script - Main (Old)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% -------------------
clear; close all; clc; format short;
%Inputs
c = .999;
sigt = 1.0;
sigs = c*sigt;
siga = sigt-sigs;
D = 1/(3*sigt);
% Quad info
quad = 4;
[qx, qw] = lgwt(quad,-1,1); nqx = length(qx);
% geometric info
xmin = 1e0; xmax = 1e0; xnum = 1e0;
ddx = get_logarithmic_x(xmin, xmax, xnum)'; nx = length(ddx);
% Phase info
np = 1;
% allocate space
max_SI    = zeros(nx, 1);
max_SIDSA = zeros(nx, 1);
SI_eigen = zeros(np,nx);
SIDSA_eigen = zeros(np,nx);

disp('-> Computing Eigen Spectrums.'); rev_str = [];
for mm=1:nx
    msg = sprintf('      -> Computing spectrum for Mesh %d of %d.',mm,nx);
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
    
    dx = ddx(mm);
    phase = linspace(0,0,np)';
%     phase = linspace(0,2*pi/dx,np)';
    e_phase = exp(1i*dx*phase);
    % build local matrices
    m = [2,1;1,2]*dx/6;
    g = [-1,-1;1,1]/2;
    k = [1,-1;-1,1]/dx;
    eg_L = [-1,0;1,0]/dx;
    eg_R = [0,-1;0,1]/dx;
    % build full matrices
    ntot = 2; zn = zeros(ntot);
    I = eye(ntot);
    T = cell(np, 1);
    S = sigs*m;
    A = D*k + siga*m;
    LL = cell(np, quad);
    cPM = cell(np, 1);
    for p=1:np
        T{p} = zn;
        cPM{p} = [1,0;0,e_phase(p)];
    end
    SS = S/2;
    % Loop through quadrature
    for q=1:nqx
        % build angle matrices
        L = sigt*m - qx(q)*g;
        for p=1:np
            PM = cPM{p};
            LL{p,q} = L*PM;
            % Based on direction
            if qx(q) > 0
                LL{p,q} = LL{p,q} + qx(q)*[0,-1;0,e_phase(p)];
            else
                LL{p,q} = LL{p,q} + qx(q)*[-1,0;e_phase(p),0];
            end
            SSPM = SS*PM;
            T{p} = T{p} + qw(q)*(LL{p,q}\SSPM);
        end
    end
    % loop through phases and calculate diffusion terms
    pp = 1; c = 2;
    kpen = max(c*pp*(pp+1)*D/dx,.25);
    for p=1:np
        PM = cPM{p};
        AA = A*PM; BB = S*PM;
        % Diffusion Boundary Terms
        % ------------------------
        e0 = 1;
        e  = e_phase(p);
        em = exp(-1i*phase(p)*dx);
        ep = exp(2i*phase(p)*dx);
        de  = [-1;1]/dx; det = [-1,1]/dx;
        bL = [1,0]; bR = [0,1];
        %Build Here...
        nL = -1; nR = 1;
        % ( [[u]] , [[b]] )
        AA = AA + kpen*[1,-1;-e,e];
        % ( {{Du}} , [[b]] )
        AA = AA + nR*D/2*[-(e+e0);(ep+e)]/dx*[0,-1];   % right boundary
        AA = AA + nL*D/2*[-(e0+em);(e0+e)]/dx*[-1,0];  % left boundary
        % ( [[u]] , {{Db}} )
        AA = AA + nL*D/2*[1/dx;-1/dx]*[e0,-e0];        % left boundary
        AA = AA + nR*D/2*[-1/dx;1/dx]*[e,-e];          % right boundary
        
        % Left Boundary
        % -------------
%         n = -1;
%         % Mass Terms
%         AA = AA + kpen*[1,-1;-1,1];
%         % Gradient Terms
%         % +,+
%         AA = AA + .5*D*n*(bR')*(det.*em);
%         AA = AA + .5*D*n*(de) *(bR.*em);
%         % -,-
%         AA = AA - .5*D*n*(bL')*(det.*epin);
%         AA = AA - .5*D*n*(de) *(bL.*epin);
%         % +,-
%         AA = AA + .5*D*n*(bR')*(det.*epin);
%         AA = AA - .5*D*n*(de) *(bL.*epin);
%         % -,+
%         AA = AA - .5*D*n*(bL')*(det.*em);
%         AA = AA + .5*D*n*(de) *(bR.*em);
%         % Right boundary
%         % --------------
%         n = 1;
%         % Mass Terms
%         AA = AA + kpen*[1,-1;-e,e];
%         % Gradient Terms
%         % +,+
%         AA = AA + .5*D*n*(bL')*(det.*ep);
%         AA = AA + .5*D*n*(de) *(bL.*ep);
%         % -,-
%         AA = AA - .5*D*n*(bR')*(det.*epin);
%         AA = AA - .5*D*n*(de) *(bR.*epin);
%         % +,-
%         AA = AA + .5*D*n*(bL')*(det.*epin);
%         AA = AA - .5*D*n*(de) *(bR.*epin);
%         % -,+
%         AA = AA - .5*D*n*(bR')*(det.*ep);
%         AA = AA + .5*D*n*(de) *(bL.*ep);
        % Calculate All Eigenvalues
        % -------------------------
        % SI
        % --
        mat_SI = T{p};
        t_SI = max(abs(eig(mat_SI)));
        SI_eigen(p,mm) = t_SI;
        if t_SI > max_SI(mm), max_SI(mm) = t_SI; end
        % SI + DSA
        % --------
        mat_SIDSA = T{p} + (AA\BB)*(T{p}-I);
        t_SIDSA = max(abs(eig(mat_SIDSA)));
        SIDSA_eigen(p,mm) = t_SIDSA;
        if t_SIDSA > max_SIDSA(mm), max_SIDSA(mm) = t_SIDSA; end
    end

end
fprintf(rev_str); rev_str = [];

% Normalize to c
max_SI_norm = max_SI/c;
max_SIDSA_norm = max_SIDSA/c;


