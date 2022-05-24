function [noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1]=fix_iteration_pro(level,noise,M,N,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy)
[nr,nc]=size(noise);
row_odd_idx(row_odd_idx<1) = 1;
row_odd_idx(row_odd_idx>nr) = nr;
col_odd_idx(col_odd_idx<1)=1;
col_odd_idx(col_odd_idx>nc)=nc;
z1=noise(row_odd_idx,col_odd_idx);
noise_ch1=zeros(M*N,size(row_odd,2)*size(col_odd,2));
% phi_idx1=zeros(M*N,size(row_odd,2)*size(col_odd,2));%reshape patch as a vertical vector
gradidx_ch1=zeros((M-1)*(N-1),size(row_odd,2)*size(col_odd,2));%record the grad in the inner patch 

nr=size(row_odd,2);
nc=size(col_odd,2);
for ii = 1:nr
    noise_ch1(:, (ii-1)*nc+1:ii*nc)= reshape(z1((ii-1)*M+1:ii*M,:), M*N, []);
%    phi_idx1(:, (ii-1)*nc+1:ii*nc)=reshape(idx_tmpw((ii-1)*M+1:ii*M,:), M*N, []);
%     gradidx_ch1(:,(ii-1)*nc+1:ii*nc)=reshape(grad_idx((ii-1)*(M-1)+1:ii*(M-1),:),(M-1)*(N-1),[]);
end
phi=ker(level);
phitmp=phi(:);
nphi_1=phitmp*ones(1,size(noise_ch1,2));
%nphi_1=nphi_1.*phi_idx1;
s1=sum(nphi_1.*nphi_1,1);
vx_1=nphi_1(bnd_idx(:,1),:)-nphi_1(bnd_idx(:,2),:);
vy_1=nphi_1(bnd_idy(:,1),:)-nphi_1(bnd_idy(:,2),:);

dv1=zeros(M*size(row_odd,2),N*size(col_odd,2));
for ii=1:size(row_odd,2)
    tmp=nphi_1(:,(ii-1)*size(col_odd,2)+1:ii*size(col_odd,2));
    dv1((ii-1)*M+1:ii*N,:)=reshape(tmp,M,[]);  
end


end