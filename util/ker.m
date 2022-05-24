
function [k_a ] =ker( level )

if level==1
    k_a=zeros(2^(level-1)+2);
else
    k_a=zeros(2^(level-1)+2);
    [m,n]=size(k_a);
    k_a=zeros(m+1,n+1);
end
k_a(2:end-1,2:end-1)=1;
end


