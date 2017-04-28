function out=ged(v)
% lambda defined by Nelson (2.5)   
lambda=sqrt(2^(-2/v)*gamma(1/v)/gamma(3/v));
% exponential ditribution parameter
theta=0.5/lambda^v;
const=exp(theta*(v^(v\(1-v))-v^(1\(1-v))));

x=-1;
r=rand;
while(x<r)
   e=-log(rand)/theta; %draw from exp
   x=const*exp(-theta*(e^v-e)); 
end

if rand>0.5
    e=-e; %symmetry
end

out=e;

end