    
% model: J(U)=\int(u-f)^2+\lambda|\nabla (u+beta)|
% f noisy image
% %u0 clean image
% % %w  result
close all
clear all

addpath('images');
addpath('util');

var=[20];%:noise level
alpha=[20];%regularization parameter
max_level=4;
load images/lena.mat
u0=im;
f=u0+randn(size(u0))*var;
 
 [ w ,Energy,Energy_out,error,error_out,t] = MMC_code(f,alpha,max_level);
psnr_u=psnr(uint8(w),uint8(u0))
ssim_u=ssim(uint8(w),uint8(u0));
figure;imshow(w,[])
