a=0.2;
M=4^9;
Nt=1000;
diff=zeros([1,Nt]);
for i=1:Nt
    xj=randp(M);
    Gxj=g(xj,a);
    diff(i)=mean(Gxj);
end
diff=diff-Iexact(a);
histogram(diff)
title('$M=4^9$ with $N_{trial}=1000$','interpreter','latex')
ylabel('Counts')
xlabel('$\hat{I}-I_{exact}$','interpreter','latex')
set(gca,'YScale','log')