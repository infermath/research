function s=montecarlo(N,T)

par =[-9.0582,(0.2163e-03)^2,0.2176,1.9132,0.9800,0.3428,0.5988;-9.1024,(0.002e-03)^2,0.1277,1.3495,0.9800,1.1317,-0.1621];
c=0.7349;

n=2;
cor=ones(n)*c; % generate correlation matrix
for i=1:n
    cor(i,i)=1;
end

omega=par(:,1);  theta=-sqrt(par(:,2));   delta= par(:,3);  thickness=par(:,4);
ma=par(:,5);  ar=par(:,6:7);

q=size(ma,2);
p=size(ar,2);
s=struct;
l=ones(1,length(T));

for i=1:N
    i
    returns=simEGARCH(cor,T(end),omega,theta,delta,thickness,ma,ar); 
    for j=1:length(T)
        [table,corr,z]=corrmatrix2(returns(1:T(j),:),q,p);
        if (unitrootcheck(table(:,end-p+1:end))==0)
            s(j).correlation(l(j),:)=corr(1,2);
            s(j).varz1(l(j),:)=corr(1,1);
            s(j).varz2(l(j),:)=corr(2,2);
            s(j).shocks=z;
            for k=1:n
                s(j).asset(k).estimator(l(j),:)=table(k,:); 
                s(j).asset(k).average=mean(s(j).asset(k).estimator);
                s(j).asset(k).stdev=std(s(j).asset(k).estimator);
            end
            l(j)=l(j)+1;
        end
    end
    save mciteration
end

for j=1:length(T)
    s(j).SampleSize=T(j);
    s(j).AverageCorrelation=nanmean(s(j).correlation);
    s(j).Averagevarz1=nanmean(s(j).varz1);
    s(j).Averagevarz2=nanmean(s(j).varz2);
    s(j).StandardDeviationOfCorrelation=nanstd(s(j).correlation);
    s(j).SkewnessCorrelation=skewness(s(j).correlation);
    s(j).KurtosisCorrelation=kurtosis(s(j).correlation);
    s(j).StandardizedCorrelation=sort(sqrt(T(j))*(s(j).correlation-c));
    s(j).AverageStandardizedCorrelation=nanmean(s(j).StandardizedCorrelation);
    s(j).StandardDeviationOfStandardizedCorrelation=nanstd(s(j).StandardizedCorrelation);
end
end