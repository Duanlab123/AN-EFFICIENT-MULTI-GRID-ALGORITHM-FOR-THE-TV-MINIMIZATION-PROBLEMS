function w=MMC_ROF(w,opt,belta,level,iter_fix)
alpha = opt.alpha;
row_odd_idx = opt.row_odd_idx;
row_even_idx = opt.row_even_idx;
col_odd_idx = opt.col_odd_idx;
col_even_idx = opt.col_even_idx;
row_odd = opt.row_odd;
row_even = opt.row_even;
col_odd = opt.col_odd;
col_even = opt.col_even;
bnd_idx = opt.bnd_idx;
bnd_idy = opt.bnd_idy;
noise=opt.noise;
s1=opt.s1;
s2=opt.s2;
s3=opt.s3;
s4=opt.s4;
noise_ch1=opt.noise_ch1;nphi_1=opt.nphi_1;vx_1=opt.vx_1;vy_1=opt.vy_1;dv1=opt.dv1;
noise_ch2=opt.noise_ch2;nphi_2=opt.nphi_2;vx_2=opt.vx_2;vy_2=opt.vy_2;dv2=opt.dv2;
noise_ch3=opt.noise_ch3;nphi_3=opt.nphi_3;vx_3=opt.vx_3;vy_3=opt.vy_3;dv3=opt.dv3;
noise_ch4=opt.noise_ch4;nphi_4=opt.nphi_4;vx_4=opt.vx_4;vy_4=opt.vy_4;dv4=opt.dv4;
gradidx_ch1=opt.gradidx_ch1;
gradidx_ch2=opt.gradidx_ch2;
gradidx_ch3=opt.gradidx_ch3;
gradidx_ch4=opt.gradidx_ch4;
[w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,belta,alpha,iter_fix,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1);
[w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_even_idx,row_odd,col_even,bnd_idx,bnd_idy,level,belta,alpha,iter_fix,noise_ch2,nphi_2,vx_2,vy_2,dv2,s2,gradidx_ch2);
[w]=MMC_DDM_ROF(noise,w,row_even_idx,col_odd_idx,row_even,col_odd,bnd_idx,bnd_idy,level,belta,alpha,iter_fix,noise_ch3,nphi_3,vx_3,vy_3,dv3,s3,gradidx_ch3);
[w]=MMC_DDM_ROF(noise,w,row_even_idx,col_even_idx,row_even,col_even,bnd_idx,bnd_idy,level,belta,alpha,iter_fix,noise_ch4,nphi_4,vx_4,vy_4,dv4,s4,gradidx_ch4);
end