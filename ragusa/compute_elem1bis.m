function [mtot,msca,g,e,NPAR]=compute_elem1bis(tot,sca,ncells)

% mass
m = [2 1; 1 2]/3;
% m=eye(2);

% jac = dx(1:ncells)/2;
mtot = tot*m;
msca = sca*m;
% gradient
g = [1 1 ; -1 -1]/2;
gr = kron( eye(ncells), g );
% edge
% if(ncells==1)
%     ed  =diag([-1;1]);
%     e{1}=ed;
%     e{1}(1,2)=e{1}(1,1);
%     e{1}(1,1)=0;
%     e{2}=ed;
%     e{2}(2,1)=e{2}(2,2);
%     e{2}(2,2)=0;   
% else
%     e{1}=[0   0  0 -1;...
%           0  +1  0  0;...
%           0  -1  0  0;...
%           0   0  0 +1];
%     e{2}=[-1  0  0  0;...
%           0   0 +1  0;...
%           0   0 -1  0;...
%          +1   0  0  0];
% end
e{1} = zeros(ncells*2,ncells*2);
e{2} = zeros(ncells*2,ncells*2);
for iel=1:ncells
    istart = (iel-1)*2+1;
    iel0 = iel-1;
    if (iel0==0) iel0 = ncells; end
    istart0 = (iel0-1)*2+1;
    iel2 = iel+1;
    if (iel2>ncells) iel2 = 1; end
    istart2 = (iel2-1)*2+1;
    e{1}(istart  ,istart0+1) = -1;
    e{1}(istart+1,istart+1 ) = 1;
    e{2}(istart  ,istart   ) = -1;
    e{2}(istart+1,istart2  ) = 1;
end
% stiffness
NPAR.k=[ 1 -1;-1  1] /2;
% DG diffusion
NPAR.intint{1} = [ 0 0 ; -1 1]/2;
NPAR.intint{2} = NPAR.intint{1}'; 
NPAR.extint{1} = [-1 1 ;  0 0]/2;
NPAR.extint{2} = [0 -1 ;  0 1]/2;
NPAR.intext{1} = [ 0 0 ; -1 1]/2;
NPAR.intext{2} = [-1 0 ;  1 0]/2;
NPAR.extext{1} = [-1 1 ;  0 0]/2;
NPAR.extext{2} = NPAR.extext{1}'; 

NPAR.cdif   = 1./(3*tot);
% 
NPAR.c=1;
NPAR.p=1;

NPAR.dcf    = [ 1 -1; -1 1];
NPAR.dcf_ii = [ 1 -1; -1 1];
NPAR.dcf_ee = [ 1 -1; -1 1];
NPAR.dcf_ie =-[ 1 -1; -1 1];
NPAR.dcf_ei =-[ 1 -1; -1 1];
    