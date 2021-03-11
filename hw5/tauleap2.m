N=10000; drs=[0.001 0.003]; d=0.003; g=0.05; a0=0.5; a1=0.01; 
T=15000;  
tau=0.1;

Xvec=NaN([3,T/tau+2]);
tvec=NaN([1,T/tau+2]);

S=zeros([6,15]); 
for i=1:3
    S(i,i)=1;
    S(i,3+i)=-1;
    S(i,6+i)=-1; S(3+i,6+i)=1;
    S(i,9+i)=1; S(3+i,9+i)=-1;
    S(3+i,12+i)=-1;
end

for i=1:2
    dr=drs(i);
    X=zeros([6,1]); t=0;
    counter=1;
    while t<T
        tvec(counter)=t;
        Xvec(:,counter)=X(1:3);
        counter=counter+1;

        rk=reaction(X,N,dr,d,g,a0,a1);
        poiz=poissrnd(tau*rk);

        Xmid=X+0.5*tau*S*poiz.';
        rk=reaction(Xmid,N,dr,d,g,a0,a1);
        poiz=poissrnd(tau*rk);
        X=X+S*poiz.';

        if any(X<0)
            disp('stop')
            break
        end
        t=t+tau;
    end

    tvec=tvec(1:counter-1);
    Xvec=Xvec(:,1:counter-1);
    
    subplot(2,1,i)
    hold on;
    plot(tvec,Xvec(1,:),'DisplayName','%X_1%')
    title(sprintf( '$d_r=%.3f$',dr));
    ylabel('Population')
    xlabel('Time')
    hold off;
end