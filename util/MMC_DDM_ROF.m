function [w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,belta,alpha,iter_fix,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1)
% [nr,nc]=size(noise);

w1=w(row_odd_idx,col_odd_idx);
%z1=noise(row_odd_idx,col_odd_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%% we have two solver[fix iteration and formula
%%%%%%%%%%%%%%%%%%%%%%%%%% solution]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dw1=MMC_fix(w1,bnd_idx,bnd_idy,row_odd,col_odd,level,belta,alpha,iter_fix,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,gradidx_ch1);
if row_odd_idx(1)==row_odd_idx(2)
dw1(1,:)=dw1(2,:);
end
if row_odd_idx(end)==row_odd_idx(end-1)
dw1(end,:)=dw1(end-1,:);
end
if col_odd_idx(1)==col_odd_idx(2)
dw1(:,1)=dw1(:,2);
end
if col_odd_idx(end)==col_odd_idx(end-1)
dw1(:,end)=dw1(:,end-1);
end
w(row_odd_idx,col_odd_idx) = w1+dw1;
end