a=0.2;
N=3:1:9;
Nt=1000;
em=zeros(size(N));
sigm=zeros(size(N));
chebm=zeros(size(N));
cltm=zeros(size(N));
for k=1:length(N)
    M=4^N(k);
    diff=zeros([1,Nt]);
    sighat=zeros([1,Nt]);
    for i=1:Nt
        xj=randp(M);
        Gxj=g(xj,a);
        diff(i)=mean(Gxj);
        sighat(i)=sqrt(var(Gxj)/M);
    end
    diff=diff-Iexact(a);
    em(k)=mean(abs(diff));
    sigm(k)=mean(sighat);
    cltm(k)=clterr(a,0.99,M);
    chebm(k)=cheberr(a,0.99,M);
end
hold on;
plot(N,em,'DisplayName','$\varepsilon_m$')
plot(N,sigm,'DisplayName','$\hat{\sigma}_{\hat{I}}$')
plot(N,cltm,'DisplayName','CLT Bound')
plot(N,chebm,'DisplayName','Chebyshev Bound')
title('$M=4^n$ and $N_{trial}=1000$','interpreter','latex')
xlabel('$4^n$','interpreter','latex')
set(gca,'YScale','log')
hl = legend('show');
set(hl, 'Interpreter','latex')
hold off;