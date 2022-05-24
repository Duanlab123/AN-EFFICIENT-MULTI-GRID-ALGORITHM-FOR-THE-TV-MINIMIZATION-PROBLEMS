function dU1=MMC_fix(w1,bnd_idx,bnd_idy,row_odd,col_even,level,belta,alpha,iter_fix,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1)
[pr,pc]=size(ker(level));
itmax=10;tauinner=1e-1;
U1ch=zeros(pr*pc,size(row_odd,2)*size(col_even,2));
%noisech=U1ch;

nr=size(row_odd,2);
nc=size(col_even,2);

for ii = 1:nr
    U1ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
 %   noisech(:, (ii-1)*nc+1:ii*nc)= reshape(z1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
end
noisech=noise_ch1;
zbar=noisech-U1ch;
nphi=nphi_1;
% phi=ker(level);
% phitmp=phi(:);
% nphi=phitmp*ones(1,size(U1ch,2));
zstar=zbar.*(nphi);
%s1=sum(nphi.*nphi,1);
zstar=sum(zstar,1);

ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);
% vx=nphi(bnd_idx(:,1),:)-nphi(bnd_idx(:,2),:);
% vy=nphi(bnd_idy(:,1),:)-nphi(bnd_idy(:,2),:);
vx=vx_1;
vy=vy_1;
c=zeros(1,size(ux,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%an-isotropic iteration%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  A1=vx.^2+vy.^2;
  a=ones(size(ux,1),1);
   C1=ux.*vx+uy.*vy;
for iter=1:itmax
    cold=c;
    tmpc=a*c;
    ucx=ux+tmpc.*vx;
    ucy=uy+tmpc.*vy;
    %          A1=A1.*tmpc;%gradient descend
    B1=ucx.^2+ucy.^2+belta;
    B1=sqrt(B1);
    A=A1./B1;
    firs=sum(A,1);
   
    C=C1./B1;
    sec=sum(C,1);
    c=(zstar-alpha.*sec)./(alpha.*firs+s1);

    if norm(cold-c,'inf')<tauinner
        break
    end
end
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));

dU1=dc1.*dv1;
end

