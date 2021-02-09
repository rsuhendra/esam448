nbins=51;
t=-1:2/nbins:1;
pj=zeros(1,nbins);
for i=1:nbins
    pj(i)=bigP(t(i+1),0.1,0.1)-bigP(t(i),0.1,0.1);
end
y=zeros(1,7);
for i=2:8
    data=reject(0.1,0.1,10^i);
    [Nj,edges] = histcounts(data,t);
    y(i-1)=chi2(pj,Nj);
end
x=2:1:8;

ylabel('\chi 2 statistic')
xlabel('Powers of 10')

plot(x,y)