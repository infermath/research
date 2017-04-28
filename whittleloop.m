function out=whittleloop(data,tol,ma,ar)
%considers GED case
%EGARCH parameters from Nelson (1991), 2 - assumed normal ditribution, ARMA parametes =0 white noise 
par=zeros(1,3+ar+ma);
par(1)=0.0139;
par(2)=0.1559;
par(3)=2;
prev=zeros(1,3+ar+ma);

    while (norm(par-prev,2)>tol)
        prev=par;
        par=whittle(data,par,ma,ar);
    end

    % omega estimation
    v=par(3);
    omega=mean(log(data.^2))-2/v*(psi(0,1/v)-log((gamma(3/v)/gamma(1/v))^(v/2)));
    out=[omega,par];
 
end