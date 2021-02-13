a=0.8;
M=10;
Nt=10000;
diff=zeros([1,Nt]);
vr=zeros([1,Nt]);
for i=1:Nt
    xj=randp(M);
    Gxj=g(xj,a);
    diff(i)=mean(Gxj);
    vr(i)=var(Gxj);
end
diff=diff-Iexact(a);
diff=abs(diff);
v=mean(vr);
cheb=cheberr(a,0.99,M,v);
clt=clterr(a,0.99,M,v);
chebperc=sum(diff<cheb)/100;
cltperc=sum(diff<clt)/100;


cltmsg=sprintf('%.2f < CLT bound',cltperc);
chebmsg=sprintf('%.2f < Chebyshev bound',chebperc);
histogram(diff)
xline(clt,'DisplayName',cltmsg)
xline(cheb,'DisplayName',chebmsg)
title('$\alpha=0.8, M=10$ with $N_{trial}=10000$','interpreter','latex')
ylabel('Counts')
xlabel('$|\hat{I}-I_{exact}|$','interpreter','latex')
legend