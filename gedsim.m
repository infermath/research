function x = gedsim(v,cor,T)

%Choleski decomposition
U=chol(cor);
n=length(v);
x=zeros(T,n);

for i=1:T
    for j=1:n
       z(j)=ged(v(j)); 
    end
    x(i,:)=z*U;
end

end

