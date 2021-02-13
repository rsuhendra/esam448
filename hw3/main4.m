a=0.2;
M=10^4;
Nt=1000;
diff=zeros([1,Nt]);
for i=1:Nt
    xj=randp(M);
    Gxj=g(xj,a);
    diff(i)=mean(Gxj);
end
diff=diff-Iexact(a);
diff=abs(diff);
cheb=cheberr(a,0.99,M);
clt=clterr(a,0.99,M);
chebperc=sum(diff<cheb)/10;
cltperc=sum(diff<clt)/10;


cltmsg=sprintf('%.1f < CLT bound',cltperc);
chebmsg=sprintf('%.1f < Chebyshev bound',chebperc);
histogram(diff)
xline(clt,'DisplayName',cltmsg)
xline(cheb,'DisplayName',chebmsg)
title('$M=10^5$ with $N_{trial}=1000$','interpreter','latex')
ylabel('Counts')
xlabel('$|\hat{I}-I_{exact}|$','interpreter','latex')
legend