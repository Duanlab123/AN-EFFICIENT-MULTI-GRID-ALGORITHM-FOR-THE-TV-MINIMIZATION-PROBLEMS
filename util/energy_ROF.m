function J_U=energy_ROF(U,I,alpha)
[m,n]=size(U);
DxU = [U(:,1:end-1) - U(:,2:end) zeros(m,1)];
DyU = [U(1:end-1,:) - U(2:end,:); zeros(1,n)];
J_U=sum(sum((U-I).^2))+alpha*sum(sum(sqrt(DxU.^2+DyU.^2)));

