function estimator = whittle(data,initial,ma,ar)

% GED parameters
mi=@(theta) gamma(2/theta(3))/(gamma(3/theta(3))*gamma(1/theta(3)))^0.5;
alpha0=@(theta) (2/theta(3))^2*psi(1,1/theta(3));
beta0=@(theta) theta(1)^2+theta(2)^2*(1-mi(theta)^2);
gamma0=@(theta) 2*theta(2)/theta(3)*mi(theta)*(psi(0,2/theta(3))-psi(0,1/theta(3)));

% transfer function phi
a=@(z,theta) 1;
for j=1:ma
    a=@(z,theta) a(z,theta)+theta(3+j)*z.^j;
end
b=@(z,theta) 1;
for j=1:ar
    b=@(z,theta) b(z,theta)-theta(3+ma+j)*z.^j;
end
phi =@(z,theta) a(z,theta)./b(z,theta);

% periodogram (T-1x1 vector)
logx2=log(data.^2);
T=length(data);
t=1:T-1;
freq=2*pi*t'/T;
I = (abs(exp(1i*freq*[t,T])*logx2).^2)/(2*pi*T);

% spectral density function f (T-1x1 vector function)
f =@(theta) 1/2/pi*(alpha0(theta)+beta0(theta)*abs(phi(exp(1i*freq),theta)).^2+ gamma0(theta)*(exp(1i*freq).*phi(exp(1i*freq),theta)+exp(-1i*freq).*phi(exp(-1i*freq),theta)));

% discrete Whittle function
Q =@(theta) sum(log(f(theta))+I./f(theta))/T;

% minimization
options = optimset('MaxFunEvals',50000,'MaxIter',50000);
epsilon=0.01;
estimator=fminsearchbnd(@(theta) Q(theta),initial,[0,0,0,-ma*ones(1,ma)+epsilon,-ar*ones(1,ar)+epsilon],[1,1,100,ma*ones(1,ma)-epsilon,ar*ones(1,ar)]-epsilon,options);

end