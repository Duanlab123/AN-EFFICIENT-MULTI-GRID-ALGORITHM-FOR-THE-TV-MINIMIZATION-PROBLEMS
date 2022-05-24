function dU1=NewDMtest(z1,w1,bnd_idx,bnd_idy,row_odd,col_even,level,belta,alpha,idxtmpw)
[pr,pc]=size(ker(level));
U1ch=zeros(pr*pc,size(row_odd,2)*size(col_even,2));
noisech=U1ch;
phiidx=U1ch;
gradidxch=zeros((pr-1)*(pc-1),size(row_odd,2)*size(col_even,2));
idxtmpw=double(idxtmpw);
%we need the domain the gradient
ctrx=2:pr;
row=zeros(1,(pr-1)*size(row_odd,2));
row(1,1:size(ctrx,2))=ctrx;
ctry=ctrx;
col=zeros(1,(pc-1)*size(col_even,2));
col(1,1:size(ctry,2))=ctry;
for i=2:size(row_odd,2)
    row(1,(i-1)*size(ctrx,2)+1:i*size(ctrx,2))=row(1,(i-2)*size(ctrx,2)+1:(i-1)*size(ctrx,2))+pr;
end
for i=2:size(col_even,2)
   col(1,(i-1)*size(ctry,2)+1:i*size(ctry,2))=col(1,(i-2)*size(ctry,2)+1:(i-1)*size(ctry,2))+pc;   
end
row=sort(row);
col=sort(col);
gradidx=idxtmpw(row,col);
nr=size(row_odd,2);
nc=size(col_even,2);
for ii = 1:nr
    U1ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    noisech(:, (ii-1)*nc+1:ii*nc)= reshape(z1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    phiidx(:, (ii-1)*nc+1:ii*nc)=reshape(idxtmpw((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    gradidxch(:,(ii-1)*nc+1:ii*nc)=reshape(gradidx((ii-1)*(pr-1)+1:ii*(pr-1),:),(pr-1)*(pc-1),[]);
end
zbar=noisech-U1ch;
phi=ker(level);
phitmp=phi(:);
nphi=phitmp*ones(1,size(U1ch,2));

zstar=zbar.*(nphi);
s1=sum(nphi.*nphi,1);
zstar=sum(zstar,1)./s1;

% vy(tmpidy1(idx2),:)=phitmp(tmpidy1(idx2),:);
ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);
vx=nphi(bnd_idx(:,1),:)-nphi(bnd_idx(:,2),:);
vy=nphi(bnd_idy(:,1),:)-nphi(bnd_idy(:,2),:);
% ux=ux.*gradidxch;
% uy=uy.*gradidxch;
% vx=vx.*gradidxch;
% vy=vy.*gradidxch;%get the gradient in the matrix

idxzeros= vx~=0;
idyzeros= vy~=0;
vx=vx(idxzeros(:,1),:);
vy=vy(idyzeros(:,1),:);
ux=ux(idxzeros(:,1),:);
uy=uy(idyzeros(:,1),:);
%get the gradient in the matrix

gradu=[ux;uy];
tmpa=-gradu./[vx;vy];
[u,idx]=sort(tmpa,1);
% u=sort(gradu,1);
tmpwi=alpha*[abs(vx);abs(vy)];
wi=zeros(size(tmpwi));
 for i=1:size(tmpwi,2)
     tmp=tmpwi(:,i);
     wi(:,i)=tmp(idx(:,i));
 end
 W_i=zeros(size(wi,1)+1,size(wi,2));
 W_i(1,:)=sum(wi,1);
 for tmpi=2:size(wi,1)+1
   W_i(tmpi,:)=W_i(tmpi-1,:)-2*wi(tmpi-1,:);  
 % W_i(tmpi)=-sum(wi(1:tmpi-1))+sum(wi(tmpi:end));
 end
 p_itmp=sign(W_i).*abs(W_i)./(ones(size(W_i,1),1)*s1);%error1
 p_i=ones(size(p_itmp,1),1)*zstar+p_itmp;

 data=[u;p_i];
 c=median(data,1);

dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
dv1=zeros(pr*size(row_odd,2),pc*size(col_even,2));
for ii=1:size(row_odd,2)
    tmp=nphi(:,(ii-1)*size(col_even,2)+1:ii*size(col_even,2));
    dv1((ii-1)*pr+1:ii*pr,:)=reshape(tmp,pr,[]);  
end
dU1=dc1.*dv1;

end