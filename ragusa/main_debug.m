clear; % clears all variables from memory
close all; % closes all plotting figures
clc; format short;

% load the angular quadrature
SNQ.sn=4;
[SNQ]=loadquadrature(SNQ);

% material
tot = 1;
sca = 0.999;

ncells=1;

% compute elementary matrix
% mt = mass matrix for total xs
% ms = mass matrix for scattering xs
% g  = gradient matrix
% e =  edge matrix e{1} for mu>0, e{2} for mu<0
% NPAR = structure containing the IP matrices
[mt_nojac,ms_nojac,g,e,NPAR]=compute_elem1bis(tot,sca,ncells);

n=40;
dx_list = 2.^((-n:n)/4);
for i=1:length(dx_list)
    disp(i);
    
    dx = dx_list(i);  
    % max lambda for periodic plot
    lambda_max=4*pi*1/dx;
    
    % compute elementary matrix
    jac = dx/2;
    mt=mt_nojac*jac;
    ms=ms_nojac*jac;
    
    % nbr of points in plot, ie, number of different lambda
    nlambda=200;
    ind=0;
    for lambda=0 : lambda_max/nlambda : lambda_max
        % T = transport matrix for all directions
        % S = transport scattering matrix
        % F = sum_{m} w_m T_m^{-1} S : matrix to obtain the scalar flux after 1 SI
        % P = F +A^{-1}(F_Id) : matrix to obtain the scalar flux with SI+DSA
        % A = DG IP matrix
        [T,S,F,P,A]=compute_T1(SNQ,lambda,mt,ms,g,e,dx,ncells,NPAR);
        ind=ind+1;
        % the results from vp and vp2 must be the same
        vp(ind) =max(abs( eig( T\S )));
        vp2(ind)=max(abs( eig( F   )));
        % preconditioned eig value
        vp3(ind)=max(abs(eig(P)));
        %    vp_(ind)=min(abs(q));
        %    vp2_(ind)=min(abs(eig(F)));
    end
    
    % plots
    x=0 : lambda_max/nlambda : lambda_max;
    
    subplot(3,1,1)
    plot(x,vp); grid on
    axis([0 max(x) 0 1]);
    
    subplot(3,1,2)
    plot(x,vp2); grid on
    axis([0 max(x) 0 1]);
    
    subplot(3,1,3)
    plot(x,vp3); grid on
    axis([0 max(x) 0 max(vp3)]);
    drawnow
    
    np=5; % start point choice
    vp_fminsearch=0;
    for k=0:np
        [lambda,vp_] = fminsearchbnd(@(x) dsa_func(x,SNQ,mt,ms,g,e,dx,ncells,NPAR),[4*k/np*pi/dx],[0],[4*pi/dx]);
        lambda_min2(k+1)=lambda;
        if abs(vp_)>vp_fminsearch
            vp_fminsearch=abs(vp_);
            lambda_min=lambda;
        end
    end
    [max(vp) max(vp2) max(vp3) vp_fminsearch]
%     [max(vp) max(vp2) max(vp3) vp_fminsearch lambda_min]
%     lambda_min2
%     pause
    y_val_plot(i) = max(vp3);
    y2val_plot(i) = vp_fminsearch;
end

figure(4); 
semilogx(dx_list,y_val_plot); hold all;
semilogx(dx_list,y2val_plot);

% % %
% % % %%%%%%%%%%%%%%%%%%
% % % % for a given choice of lambda, get the matrices and perform
% % % % by hand SI and SI+DSA
% % % % we put a flat external source q in the volume
% % %
% % % lambda=0.0;
% % % [T,S,F,P,A]=compute_T1(SNQ,lambda,mt,ms,g,e,dx,ncells,NPAR);
% % % % FEM representation of the angular source
% % % qa=kron(dx(1:ncells), ones(1,2*SNQ.sn)) /2;qa=qa';
% % % % init
% % % psi_old=rand(2*ncells*SNQ.sn,1);
% % % % number of SI terations
% % % nbr_SI=5;
% % % % SI based on the angular flux
% % % for k=1:nbr_SI
% % %     psi_new = inv(T) * ( S*psi_old + qa/2);
% % %     err     = norm(psi_new-psi_old,2);
% % %     if(k>1)
% % %         sr = err/err_0; % estimation of the spectral radius
% % %     else
% % %         err_0=99;
% % %         sr=-1;
% % %     end
% % %     disp(sprintf('SI-v1 %i, err=%g, sr=%g',k,err,sr));
% % %     psi_old = psi_new;
% % %     err_0=err;
% % %     if(err<1e-11)
% % %         break
% % %     end
% % % end
% % % disp(' ')
% % %
% % % % FEM representation of the source
% % % q=kron(dx(1:ncells), ones(1,2)) /2; q=q';
% % %
% % % % SI based on the scalar flux
% % % % the error from both version should behave exactly the same
% % % phi_old=rand(2*ncells,1);
% % % phi_new=phi_old;
% % % for k=1:nbr_SI
% % %
% % %     phi_new = F*( phi_old +2*inv(ms)*q/2 );
% % %     err = norm(phi_new-phi_old,2);
% % %     if(k>1)
% % %         sr = err/err_0; % estimation of the spectral radius
% % %     else
% % %         err_0=99;
% % %         sr=-1;
% % %     end
% % %     disp(sprintf('SI-v2 %i, err=%g, sr=%g',k,err,sr));
% % %     phi_old=phi_new;
% % %     err_0=err;
% % %     if(err<1e-11)
% % %         break
% % %     end
% % %
% % % end
% % % disp(' ')
% % %
% % %
% % % % SI + DSA
% % % phi_old=rand(2*ncells,1)*888888888888888;
% % % err_0=99;
% % % q(1) = q(1) + 0.1;
% % % for k=1:3000
% % %
% % %     phi_1   = F*(phi_old+2*inv(ms)*q/2);
% % %     delphi  = inv(A)*ms*(phi_1-phi_old);
% % %     phi_new = phi_1 + delphi;
% % %     err = norm(phi_new-phi_old,2);
% % %     sr = err/err_0;
% % %     disp(sprintf('SI-DSA %i, err=%g, sr=%g',k,err,sr));
% % %     phi_old=phi_new;
% % %     err_0=err;
% % %     if(err<1e-11)
% % %         break
% % %     end
% % %
% % % end
% % % disp(' ')
% % %
% % % A\q
% % % abso=1./(tot-sca)
% % %
% % % max(vp)
% % % max(vp2)
% % % max(vp3)
% % %
% % % disp(sprintf('optical thickness %g',tot.*dx));
% % % % [v d]=eig(P)
% % %
% % % %
% % % % figure(2)
% % % % L=linspace(0.001,10,100);
% % % % c=1;
% % % % alpha=dot(L,L)/3+(1-c);
% % % % clear i;
% % % % for k=1:(length(L))
% % % %     num=1./(1+i*L(k)*SNQ.mu);
% % % %     om(k)=c/2*dot(SNQ.w,num);
% % % % end
% % % % e=abs(om+c/alpha*(om-1));
% % % % ee=c./L.*atan(L);
% % % % plot(L,e,L,ee,'r')
% % % % axis([0,max(L),0,1]);
