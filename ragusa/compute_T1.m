function [T,S,F,Prec,A]=compute_T1(snq,lambda,mt,ms,gg,e,dx,ncells,NPAR)
clear i;

% create the phase shift arrays
phase = [1 exp(1i*lambda*dx(1))];
% if ncells=2, increase phase size
for iel=2:ncells
    phase = [phase phase(end) phase(end)*exp(1i*lambda*dx(iel))];
end
P = diag(phase);

% compute the sum of the SN weights
sw=sum(snq.w);

% dimension of the linear system 
ndofs = ncells*2;
F=zeros(ndofs,ndofs);
TT=F;

% build the action of all transport sweeps on the scattering term
ndir = snq.sn;
for idir=1:ndir
    m=snq.mu(idir); % shortcut
    % choose the edge contribution based on the sweeping order (L to R or R to L)
    if(m>0)
        edg=e{1};
    else
        edg=e{2};
    end

    i1=(idir-1)*ndofs+1;
    i2=idir*ndofs;

    % transport matrix for tthis direction
    % = mu*gradient + total mass + edge
    % note that the jacobian dx/2 is already included in the mass matrices
    % and that it does not appear in the gradient and edge matrices for 1D 
    TT = (m*gg+mt)*P + m*P*edg;
    % one way is to build the entire SN matrix T
    % and the scattering matrix S
    T(i1:i2,i1:i2) = TT;
    S(i1:i2,:) = kron(snq.w, (ms)*P ) /sw;
    % the other way is to perform the sweep and collapse in angle to
    % obtain directly an error matrix in the scalar errors
    F = F + snq.w(idir) * inv(TT)* (ms*P/sw);
end

% DSA matrix
A=zeros(ndofs,ndofs);
jac = dx/2;
adif=NPAR.cdif ./ jac ;

% create stiff+mass matrix
% the mt ms matrices have already been multiplied by the Xs and the
% jacobian
% adif now contains the diffusion coefficient divided by the jacobian
A0 = kron( diag(adif(1:ncells)),NPAR.k ) + (mt-ms);
A=A0*P;

%%%% simpler way of getting diffusion matrix
AA=A;
% poly order
p=NPAR.p;
% constant in the penalty formula
c=NPAR.c;
ddx = NPAR.cdif(1)/dx(1);
pen = max(0.25 , c*p*(p+1)*(ddx+ddx));
% phase M--0xx1--2
pb0=P(1,1);
pb1=P(2,2);
pb2=pb1*pb1;
pbM=1/pb1;
% add penalty contribution
aux = [ 1*pb0  -1*pb0 ;...
       -1*pb1   1*pb1   ]; 
AA=AA + pen*aux;
aux = ddx/2 * [ -1*(pbM+pb0)       1*(pb0+pb1) ; ...
                 1*(pb0+pb1)      -1*(pb1+pb2) ];
AA=AA + aux;
aux = ddx/2 * [ -1*(pb0+pb1)       1*(pb0+pb1) ; ...
                 1*(pb0+pb1)      -1*(pb0+pb1) ];
AA=AA + aux;

% connectivity array
g=zeros(2,1);
gn=g;
% % %%%%%%%%%%
E111 = [0 0; 0 1];
E122 = [1 0; 0 0];
E112 = [0 0; -1 0];
E121 = [0 -1; 0 0];
E211 = 1/2*[0 0; 1 -1];
E212 = E211;
E222 = 1/2*[-1 1; 0 0];
E221 = E222;
E311 = [-1 1; 1 -1];
E322 = E311;
E312 = -[-1 1; 1 -1];
E321 = E312;
for iel=1:ncells
    istart1 = (iel-1)*2+1;
    iend1 = istart1+1;
    p1 = P(istart1:iend1,istart1:iend1);
    ddx1 = NPAR.cdif(iel)/dx(iel);
    if (iel==1)
        p0 = diag([exp(-1i*lambda*dx(end)),1]);
        ddx0 = NPAR.cdif(end)/dx(end);
        istart0 = (ncells-1)*2+1;
        iend0 = istart0+1;
        iel0 = ncells;
    else
        ddx0 = NPAR.cdif(iel-1)/dx(iel-1);
        istart0 = (iel-2)*2+1;
        iend0 = istart0+1;
        p0 = P(istart0:iend0,istart0:iend0);
        iel0 = iel-1;
    end
    if (iel==ncells)
        p2 = diag(P(end,end)*[1, exp(1i*lambda*dx(1))]);
        ddx2 = NPAR.cdif(1)/dx(1);
        istart2 = 1;
        iend2 = istart2+1;
        iel2 = 1;
    else
        ddx2 = NPAR.cdif(iel+1)/dx(iel+1);
        istart2 = (iel)*2+1;
        iend2 = istart2+1;
        p2 = P(istart2:iend2,istart2:iend2);
        iel2 = iel+1;
    end
    c1 = 0.25; c2 = 0.25;
    cc1 = 9/8; cc2 = 9/8;
     c1 = max(0.25 , c*p*(p+1)*(ddx0+ddx1));
     c2 = max(0.25 , c*p*(p+1)*(ddx2+ddx1));
     cc1 = 0; cc2 = 0;
%    c1=0; c2=0;
%    p0 = eye(2); p1=p0; p2=p1;
    A(istart1:iend1,istart1:iend1) = A(istart1:iend1,istart1:iend1) + (c1*E111+c2*E122)*p1;
    A(istart1:iend1,istart2:iend2) = A(istart1:iend1,istart2:iend2) + (c2*E112)*p2;
    A(istart1:iend1,istart0:iend0) = A(istart1:iend1,istart0:iend0) + (c1*E121)*p0;
    A(istart1:iend1,istart1:iend1) = A(istart1:iend1,istart1:iend1) + ddx1*(E211+E222)*p1;
    A(istart1:iend1,istart2:iend2) = A(istart1:iend1,istart2:iend2) + ddx2*(E212)*p2;
    A(istart1:iend1,istart0:iend0) = A(istart1:iend1,istart0:iend0) + ddx0*(E221)*p0;
    A(istart1:iend1,istart1:iend1) = A(istart1:iend1,istart1:iend1) + ddx1*(E211'+E222')*p1;
    A(istart1:iend1,istart2:iend2) = A(istart1:iend1,istart2:iend2) + ddx1*(E221')*p2;
    A(istart1:iend1,istart0:iend0) = A(istart1:iend1,istart0:iend0) + ddx1*(E212')*p0;
    A(istart1:iend1,istart1:iend1) = A(istart1:iend1,istart1:iend1) + ddx1*ddx1*(E311*cc1+E322*cc2)*p1;
    A(istart1:iend1,istart2:iend2) = A(istart1:iend1,istart2:iend2) + ddx1*ddx2*(cc2*E312)*p2;
    A(istart1:iend1,istart0:iend0) = A(istart1:iend1,istart0:iend0) + ddx1*ddx0*(cc1*E321)*p0;
end
% if(ncells==1)
% p1 = P;
% p0 = diag([ exp(-i*lambda*dx(1)), 1]);
% p2 = diag([ exp(i*lambda*dx(1)), exp(2*i*lambda*dx(1))]);
% ddx = NPAR.cdif/dx;
% c = 0.25;
% 
% % cdifp = NPAR.cdif(1)  /2;
% % cdifm = cdifp;
% % dxp = dx(1);
% % dxm = dxp;
% % c = max(0.25 , 2*c*p*(p+1)*(cdifp/dxp+cdifm/dxm) )
% 
% cc = 9/8;
% % cc=0;
% E111 = c*[0 0; 0 1];
% E122 = c*[1 0; 0 0];
% E112 = c*[0 0; -1 0];
% E121 = c*[0 -1; 0 0];
% E211 = 1/2*[0 0; 1 -1];
% E212 = E211;
% E222 = 1/2*[-1 1; 0 0];
% E221 = E222;
% E311 = cc*[-1 1; 1 -1];
% E322 = E311;
% E312 = -cc*[-1 1; 1 -1];
% E321 = E312;
% 
% A = A + (E111+E122)*p1;
% A = A + (E112)*p2;
% A = A + (E121)*p0;
% A = A + ddx*(E211+E222)*p1;
% A = A + ddx*(E212)*p2;
% A = A + ddx*(E221)*p0;
% A = A + ddx*(E211'+E222')*p1;
% A = A + ddx*(E221')*p2;
% A = A + ddx*(E212')*p0;
% A = A + ddx*ddx*(E311+E322)*p1;
% A = A + ddx*ddx*(E312)*p2;
% A = A + ddx*ddx*(E321)*p0;
% 
% else
% 
% end
% % % %%%%%%%%%%

% for iel=1:ncells
%     
%     % the 1/2 comes from D/2 in the formulae ( the averaged quantities { } in DG-IP )
%     if(ncells==1)
%         cdifp = NPAR.cdif(1)  /2;
%         cdifm = cdifp;
%         dxp = dx(1);
%         dxm = dxp;
%         g  = [1 2];
%         gn = [1 2];
%         % phase matrices
%         % we are considering a single cell, for which the phase matrix is p1;
%         % for its left element, we use p0 and for its right element, we use
%         % p2;
%         % the origin x=0 is taken at the left vertex of the element under
%         % consideration
%         p0 = diag([ exp(-i*lambda*dx(1)), 1]);
%         p1 = P;
%         p2 = diag([ exp(i*lambda*dx(1)), exp(i*lambda*(2*dx(1)))] ); 
% 
%     else
%         if(iel==1)
%             cdifp = NPAR.cdif(1)  /2;
%             cdifm = NPAR.cdif(2)  /2;
%             dxp = dx(1);
%             dxm = dx(2);
%             g  = [1 2];
%             gn = [3 4];
%             p0 = diag([ exp(-i*lambda*dx(2)), 1]);
%             p1 = P(1:2,1:2);
%             p2 = P(3:4,3:4);
%             
%         elseif(iel==2)
%             cdifp = NPAR.cdif(2)  /2;
%             cdifm = NPAR.cdif(1)  /2;
%             dxp = dx(2);
%             dxm = dx(1);
%             g  = [3 4];
%             gn = [1 2];
%             p0 = P(1:2,1:2); %diag([ exp(-i*lambda*dx(1)), 1]);
%             p2 = diag([ exp(i*lambda*(dx(1)+dx(2))), ...
%                 exp(i*lambda*(2*dx(1)+dx(2)))] );
%             p1 = P(3:4,3:4);
% 
%         else
%             error('qqq')
%         end
%     end
%     pen   = max(0.25 , 2*c*p*(p+1)*(cdifp/dxp+cdifm/dxm) );
%     pen=0.25;
%     disp(sprintf('penalty coef = %g',pen));
% 
%     
%     sg=-1;
%     % account for jacobian transformation in the grad
%     cdifp = cdifp * 2/dxp ;
%     cdifm = cdifm * 2/dxm ;
% 
%     mii = sg*( cdifp*NPAR.intint{1} + cdifp*NPAR.intint{2} ) ;
%     mei = sg*(-cdifp*NPAR.extint{1} + cdifm*NPAR.extint{2} ) ;
% 
%     mie = sg*( cdifm*NPAR.intext{1} - cdifp*NPAR.intext{2} ) ;
%     mee = sg*(-cdifm*NPAR.extext{1} - cdifm*NPAR.extext{2} ) ;
%     
%     % add the penalty in appropriate locations
%     mii( 2, 2 ) =  mii( 2, 2 ) + pen ;
%     mei( 1, 2 ) =  mei( 1, 2 ) - pen ;
% 
%     mie( 2, 1 ) =  mie( 2, 1 ) - pen ;
%     mee( 1, 1 ) =  mee( 1, 1 ) + pen ;
% 
%     % apply phase
%     mii = mii * p1 ;
%     mee = mee * p1 ;
% 
%     % left
%     mei = mei * p0 ;
%     % right
%     mie = mie * p2 ;
%     
%     A( g(:) ,g(:)  ) = A( g(:) ,g(:)  ) + mii;
%     A( gn(:),gn(:) ) = A( gn(:),gn(:) ) + mee;
% 
%     A( gn(:),g(:)  ) = A( gn(:),g(:)  ) + mei;
%     A( g(:) ,gn(:) ) = A( g(:) ,gn(:) ) + mie;
%     
%     dcf=1;
%     if(dcf)
%         taup = 1/(3*cdifp) ;
%         taum = 1/(3*cdifm) ;
%         tau  = ( taup + taum ) /2 ;
%         cc = 16;
%         pen_dcf = -min(tau/cc,9/16);
%         pen_dcf = -9/8;
% %         pen_dcf=1e9;
%         disp(sprintf('tau =%g, penalty DCF coef = %g',tau,pen_dcf));
% 
%         A( g(:) ,g(:)  ) = A( g(:) ,g(:)  ) + pen_dcf*cdifp*cdifp*NPAR.dcf*p1 ;
%         A( gn(:),gn(:) ) = A( gn(:),gn(:) ) + pen_dcf*cdifm*cdifm*NPAR.dcf*p1 ;
%         A( gn(:),g(:)  ) = A( gn(:),g(:)  ) - pen_dcf*cdifm*cdifp*NPAR.dcf*p0 ;
%         A( g(:) ,gn(:) ) = A( g(:) ,gn(:) ) - pen_dcf*cdifp*cdifm*NPAR.dcf*p2 ;
%     end
% end

% error matrix with DSA
  Prec_old = F + inv(A)*ms*P*( F- eye(ndofs) );
  Prec = F + inv(AA)*ms*P*( F- eye(ndofs) );
% Prec = inv(P)*F + inv(A)*ms*( F-eye(n)*P );

