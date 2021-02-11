a=0.2;
N=3:1:9;
N=2.^N;
Nt=1000;
em=zeros(size(N));
sigm=zeros(size(N));
for k=1:length(N)
    M=N(k);
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
end
plot(N,em,'DisplayName','em')
hold on;
plot(N,sigm,'DisplayName','sigm')
title('$M=4^3$ with $N_{trial}=1000$','interpreter','latex')
ylabel('$\varepsilon_m$','interpreter','latex')
xlabel('$M$','interpreter','latex')
set(gca,'YScale','log')
set(gca,'XScale','log')
legend;
hold off;