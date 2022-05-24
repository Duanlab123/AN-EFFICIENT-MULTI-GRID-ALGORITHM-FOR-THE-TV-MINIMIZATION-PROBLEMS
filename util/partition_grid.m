function opt=partition_grid(f,patchsize,alpha,level)
[nr,nc]=size(f);
M = patchsize(1); % number of rows in each patch
N = patchsize(2); % number of cols in each patch
pr=M-2;pc=N-2;
j=level;
npr=floor((nr-1)/2^(j-1)+1);%number of row of patch
npc=floor((nc-1)/2^(j-1)+1);%number of col of patch
    row_odd = 1:2:npr;
    row_even= 2:2:npr;
    col_odd=1:2:npc;
    col_even=2:2:npc;
row_odd_ctrx=double(1+(row_odd-1)*2^(j-1));%odd patches location
row_even_ctrx=1+(row_even-1)*2^(j-1);
col_odd_ctrx=1+(col_odd-1)*2^(j-1);
col_even_ctrx=1+(col_even-1)*2^(j-1);
pr=pr+2;pc=pc+2;


row_odd_idx = zeros(pr*length(row_odd),1);
for ii=1:length(row_odd)
row_odd_idx((ii-1)*pr+1:ii*pr,1)=row_odd_ctrx(ii)-(pr-1)/2:row_odd_ctrx(ii)+(pr-1)/2;
end
row_odd_idx(row_odd_idx<1) = 1;
row_odd_idx(row_odd_idx>nr) = nr;
row_even_idx = zeros(pr*length(row_even),1);
for ii=1:length(row_even)
    row_even_idx((ii-1)*pr+1:ii*pr,1) = row_even_ctrx(ii)-(pr-1)/2:row_even_ctrx(ii)+(pr-1)/2;
end
row_even_idx(row_even_idx<1) = 1;
row_even_idx(row_even_idx>nr) = nr;
col_odd_idx = zeros(pc*length(col_odd),1);
for ii=1:length(col_odd)
    col_odd_idx((ii-1)*pc+1:ii*pc,1) = col_odd_ctrx(ii)-(pc-1)/2:col_odd_ctrx(ii)+(pc-1)/2;
end
col_odd_idx(col_odd_idx<1) = 1;
col_odd_idx(col_odd_idx>nc) = nc;
col_even_idx = zeros(pc*length(col_even),1);
for ii=1:length(col_even)
    col_even_idx((ii-1)*pc+1:ii*pc,1) = col_even_ctrx(ii)-(pc-1)/2:col_even_ctrx(ii)+(pc-1)/2;
end
col_even_idx(col_even_idx<1) = 1;
col_even_idx(col_even_idx>nc) = nc;
pr=pr-2;pc=pc-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = pr; N = pc;
bnd_idx = zeros(4*M, 2);
bnd_idy = zeros(4*M,2);
if M==1 && N==1
    bnd_idy(:,1)=[5;8;6;9];
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)-(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)-1;
else
    bnd_idy(1:M+1,1)=(M+2)+2:(M+2):(M+2)+2+(M+2)*(N);
    bnd_idy(M+2:M+2+M,1)=bnd_idy(1:M+1,1)+N;
    bnd_idy(2*M+3:2*M+1+M,1)=(M+2)+3:(M+2)*2-1;
    bnd_idy(3*M+2:4*M,1)=bnd_idy(2*M+3:2*M+1+M,1)+(N)*(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)-1;
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)-(M+2);
end
%%%%%%%


%%%%%%%%%%%
ntau=pr*pc;
noise=f;
pr=M+2;pc=N+2;
[noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1]=fix_iteration_pro(level,noise,pr,pc,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy);
[noise_ch2,nphi_2,vx_2,vy_2,dv2,s2,gradidx_ch2]=fix_iteration_pro(level,noise,pr,pc,row_odd,col_even,row_odd_idx,col_even_idx,bnd_idx,bnd_idy);
[noise_ch3,nphi_3,vx_3,vy_3,dv3,s3,gradidx_ch3]=fix_iteration_pro(level,noise,pr,pc,row_even,col_odd,row_even_idx,col_odd_idx,bnd_idx,bnd_idy);
[noise_ch4,nphi_4,vx_4,vy_4,dv4,s4,gradidx_ch4]=fix_iteration_pro(level,noise,pr,pc,row_even,col_even,row_even_idx,col_even_idx,bnd_idx,bnd_idy);
opt.gradidx_ch1=gradidx_ch1;
opt.gradidx_ch2=gradidx_ch2;
opt.gradidx_ch3=gradidx_ch3;
opt.gradidx_ch4=gradidx_ch4;
opt.noise_ch1=noise_ch1;
opt.s1=s1;
opt.s2=s2;
opt.s3=s3;
opt.s4=s4;
opt.nphi_1=nphi_1;
opt.vx_1=vx_1;
opt.vy_1=vy_1;
opt.dv1=dv1;
opt.noise_ch2=noise_ch2;
opt.nphi_2=nphi_2;
opt.vx_2=vx_2;
opt.vy_2=vy_2;
opt.dv2=dv2;
opt.noise_ch3=noise_ch3;
opt.nphi_3=nphi_3;
opt.vx_3=vx_3;
opt.vy_3=vy_3;
opt.dv3=dv3;
opt.noise_ch4=noise_ch4;
opt.nphi_4=nphi_4;
opt.vx_4=vx_4;
opt.vy_4=vy_4;
opt.dv4=dv4;
opt.ntau=ntau;
opt.npr=npr;
opt.npc=npc;
% opt.g = g;
opt.noise=f;
opt.alpha = alpha;
opt.row_odd_idx = row_odd_idx;
opt.row_even_idx = row_even_idx;
opt.col_odd_idx = col_odd_idx;
opt.col_even_idx = col_even_idx;
opt.row_odd = row_odd;
opt.row_even = row_even;
opt.col_odd = col_odd;
opt.col_even = col_even;
opt.patchsize = patchsize;
opt.bnd_idx = bnd_idx;
opt.bnd_idy = bnd_idy;
end