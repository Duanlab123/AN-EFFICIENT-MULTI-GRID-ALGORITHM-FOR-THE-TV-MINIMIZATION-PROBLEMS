function [ w ,Energy,Energy_out,error,error_out,t] = MMC_code(z,alpha,max_level)
%parameter
w=zeros(size(z));
beta=1e-6;
iter_fix=10;
num_outer=1000;
num_iter=100;
Energy= zeros(1000*max_level,1);
J_U   = zeros (1000*max_level,1);
error = zeros(1000*max_level,1);
error_out=zeros(1000*max_level,1);
J0=energy_ROF(w,z,alpha);
J_U(1)=J0;
Energy(1)=J0;
count=2;iter=2;
%divide the grid
opt0=cell(max_level,1);
for level=1:max_level
   patchsize=size(ker(level));
   opt0{level}=partition_grid(z,patchsize,alpha,level);
end
max_loop=1;
 t1=clock;
for loop=1:max_loop
    disp(['loop number=' num2str(loop)]);  
for outer_iter=1:num_outer
    w_old_out=w;
    disp(['outer iteration' num2str(outer_iter)]);
for level=max_level:-1:1
    
    disp(['level' num2str(level)]);
    opt=opt0{level};
%     J_tmp=zeros(num_iter,1);
    error_tmp=zeros(num_iter,1);
for ii=1:num_iter  
  
    w_old=w;
    [ w ] =MMC_ROF(w,opt,beta,level,iter_fix);
    error(count-1)=norm(w-w_old,'fro')/norm(w,'fro');
    error_tmp(ii)=error(count-1);
    count=count+1;
    if error(count-2)<1e-2
          break
    end  
end
 disp(['iteration number=' num2str(ii)]);
end
J_new =energy_ROF(w,z,alpha);
J_U(iter)=J_new;iter=iter+1;
error_out(outer_iter)=norm(w-w_old_out,'fro')/norm(w,'fro');
if error_out(outer_iter)<1e-4
    break
end
end
end
   t2=clock;
   t=etime(t2,t1);
%  Energy=Energy(1:count-1);
 error=error(1:count-2);
 error_out=error_out(1:outer_iter);
Energy_out=J_U(1:outer_iter+1);
end





