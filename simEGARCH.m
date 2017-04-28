function returns=simEGARCH(cor,T,omega,theta,delta,thickness,ma,ar)

n=length(cor); % number of assets
p=size(ar,2);
q=size(ma,2);
l=max(p,q+1);

% simulating correlated schocks 
z = gedsim(thickness, cor, T+l);
err=zeros(T+l,n);
h=zeros(T+l,n);

% innovation fuction with theoretical moments
for j=1:n
    err(:,j) = theta(j)*z(:,j) + delta(j)*(abs(z(:,j))-gamma(1)/sqrt(gamma(3/2)*gamma(1/2)));
    h(:,j)=h(:,j)+omega(j); 
    if p>0
        arinv(j,:)=ar(j,end:-1:1);
    else
        arinv(j,:)=zeros(1,0);
    end
    if q>0
        mainv(j,:)=ma(j,end:-1:1);
    else
        mainv(j,:)=zeros(1,0);
    end
    const(j)=(1-sum(arinv(j,:)))*omega(j);

    for i=l+1:T+l % ARMA representation of h function   
        h(i,j)=const(j)+arinv(j,:)*h(i-p:i-1,j)+err(i-1,j)+mainv(j,:)*err(i-1-q:i-2,j);
    end
end

returns = z(l+1:end,:).*exp(0.5*h(l+1:end,:));
end