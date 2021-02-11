a=0.2;
M=10^5;
Nt=1000;
Ihat=zeros([1,Nt]);
sighat=zeros([1,Nt]);
for i=1:Nt
    xj=randp(M);
    Gxj=g(xj,a);
    Ihat(i)=mean(Gxj);
    sighat(i)=sqrt(var(Gxj)/M);
end
histogram(Ihat)
title('$M=10^5$ with $N_{trial}=1000$','interpreter','latex')
ylabel('Counts')
xlabel('$\hat{I}$','interpreter','latex')