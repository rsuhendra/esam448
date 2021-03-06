a=0.8;
M=10^5;
Nt=10000;
Ihat=zeros([1,Nt]);
sighat=zeros([1,Nt]);
for i=1:Nt
    xj=randp(M);
    Gxj=g(xj,a);
    Ihat(i)=mean(Gxj);
    sighat(i)=sqrt(var(Gxj)/M);
end
histogram(Ihat)
title('$\alpha=0.8, M=10^5$ with $N_{trial}=10000$','interpreter','latex')
ylabel('Counts')
xlabel('$\hat{I}$','interpreter','latex')

% histogram(sighat)
% title('$\alpha=0.8, M=10^5$ with $N_{trial}=10000$','interpreter','latex')
% ylabel('Counts')
% xlabel('$\hat{\sigma}_{\hat{I}}$','interpreter','latex')

set(gca,'YScale','log')

