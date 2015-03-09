function sp = dsa_func(lambda,SNQ,mt,ms,g,e,dx,ncells,NPAR)


[T,S,F,P,A]=compute_T1(SNQ,lambda,mt,ms,g,e,dx,ncells,NPAR);
% the results from vp and vp2 must be the same
vp =max(abs( eig( T\S )));
vp2=max(abs( eig( F   )));
% preconditioned eig value
vp3=max(abs(eig(P)));

if abs(vp-vp2) > sqrt(eps)
    vp
    vp2
    error('vp and vp2 are too different')
end

sp=vp3;
sp = -sp;