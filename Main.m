%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          1D LD script
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
addpath(genpath('export_fig'));
out_dir = 'outputs/';
% Inputs
plotting_bool = false;
TransportType = 'SI';
DSAType = 'MIP';
legend_position = 'northwest';
% XS inputs
c = .9;
data.IPConstant = 4;
data.XS.sigt = 1.0;
data.XS.sigs = c*data.XS.sigt;
data.XS.siga = data.XS.sigt-data.XS.sigs;
data.XS.D = 1/(3*data.XS.sigt);
% geometric info
ncells = 1; ndofs = 2*ncells;
xmin = 1e0; xmax = 1e0; xnum = 1e0;
[ddx, xtick_spacing] = get_logarithmic_x(xmin, xmax, xnum); 
ddx = ddx'; nx = length(ddx);
data.ncells = ncells; data.ndofs = ndofs;
% Quad info
quad = [2]; nquad = length(quad);
data.SN.norm = 2;
% Phase info
np = 4;
phase = linspace(0,2*pi,np)';
% Allocate Memory Space
% ---------------------
max_SI = zeros(nx,nquad);
max_SIDSA = zeros(nx, nquad);
ms_lambda = zeros(np, nx, nquad);
ms_val    = zeros(np, nx, nquad);
ms_count  = zeros(np, nx, nquad);
% Loop through Quadrature Sn Orders
% ---------------------------------
func_name = ['search_func_',TransportType,'_',DSAType];
search_func = str2func(func_name);
disp('-> Computing Eigen Spectrums.'); rev_str = [];
for q=1:nquad
    [qx, qw] = lgwt(quad(q),-1,1); nqx = length(qx);
    data.SN.num_dirs = nqx;
    data.SN.w = qw;
    data.SN.mu = qx;
    % Loop through Meshes
    for m=1:nx
        msg = sprintf('      -> Computing spectrum for Mesh %d of %d and quad %d of %d.',m,nx,q,nquad);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
        % Update dimension and phase
        data.dx = ddx(m); dx = data.dx;
        data.mats = get_1D_mats(data.dx);
        % Loop through Phases
        for p=1:np
            pp = p-1;
            [ms_lambda(p,m,q),ms_val(p,m,q),ef,out] = fminsearchbnd(@(x) search_func(x,data),2*pi/dx*(pp/np),0,2*pi/dx);
            ms_count(p,m,q) = out.funcCount;
            if abs(ms_val(p,m,q)) > max_SIDSA(m,q)
                max_SIDSA(m,q) = abs(ms_val(p,m,q)); 
            end
        end
    end
end
fprintf(rev_str); rev_str = [];

% Build output Folder if it doesn't exist
% ---------------------------------------
if ~isequal(exist(out_dir, 'dir'),7)
    mkdir(out_dir);
end
% Plot Solutions
% --------------
if plotting_bool
    fsize = 22;
    figure(1);hFig = figure(1);
    set(hFig,'Position',[1,1,1200,700])
    for q=1:nquad
        legendInfo{q} = ['S_{',num2str(quad(q)),'}'];
    end
    set(gca,'ColorOrder',[0,0,0;1,0,0;0,0,1;.5,0,.9;0,.5,0;0,1,1;1,0,1]); hold on;
    plot(ddx,max_SIDSA,'LineWidth',2.0);
    xlim([min(ddx),max(ddx)]);
    set(gca,'xscale','log');
    set(gca,'XTick',xtick_spacing);
    set(gca,'XGrid','on','XMinorGrid','off');
    set(gca,'YGrid','on','YMinorGrid','off');
    box on;
    legend(legendInfo,'Location',legend_position);
    set(gca,'FontName','Times New Roman','FontSize',fsize);
    xlabel('mfp', 'FontName', 'Times New Roman', 'FontSize', fsize, 'FontWeight', 'bold');
    ylabel('Spectral Radius', 'FontName', 'Times New Roman', 'FontSize', fsize, 'FontWeight', 'bold');
    hold off;
    cFig = figure(1);
    % Save Figures
    % ------------
    f_name = [out_dir,'DSA_1D_',TransportType,'_',DSAType,'_C=',num2str(data.IPConstant)];
    savefig(cFig,f_name)
    export_fig(f_name,'-transparent')
end