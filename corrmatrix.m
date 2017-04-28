function [par,c,z]=corrmatrix(returns,q,p)

T=size(returns,1);
n=size(returns,2);
l=max(p,q+1);
% Whittle estimation
for j=1:n
    par(j,:)=whittleloop(returns(:,j),0.1,q,p);
    arinv(j,:)=par(j,end:-1:end+1-p);
    mainv(j,:)=par(j,end-p:-1:5);
    ss(j)=sum(arinv(j,:));
end

% schocks calculation
h=zeros(T,n);
z=zeros(T,n);
err=zeros(T,n);
for i=1:T
    for j=1:n
        h(i,j)=par(j,1).*(1-ss(j));
        if (i>l)
            h(i,j)=h(i,j)+ arinv(j,:)*h(i-p:i-1,j)+err(i-1,j)+mainv(j,:)*err(i-1-q:i-2,j);
        end
    end
    z(i,:)=returns(i,:).*exp(-0.5*h(i,:)); 
    for j=1:n
        err(i,j)=par(j,2)*z(i,j)+par(j,3)*(abs(z(i,j))-gamma(1)/(gamma(1.5)*gamma(0.5))^0.5);
    end
end

% sample covariance matrix
c=(z'*z)/T;

