set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

Ns=[1 1 10 10]; drs=[0.001 0.003 0.001 0.003]; d=0.003; g=0.05; a0=0.5; a1=0.01; 
T=20000;


S=zeros([6,15]); 
for i=1:3
    S(i,i)=1;
    S(i,3+i)=-1;
    S(i,6+i)=-1; S(3+i,6+i)=1;
    S(i,9+i)=1; S(3+i,9+i)=-1;
    S(3+i,12+i)=-1;
end

for i=1:4
    dr=drs(i);
    N=Ns(i);
    X=zeros([6,1]);
    Xvec=NaN([3,30000]);
    tvec=NaN([1,30000]);
    t=0;
    counter=1;
    while t<T
        tvec(counter)=t;
        Xvec(:,counter)=X(1:3);
        counter=counter+1;

        rk=reaction(X,N,dr,d,g,a0,a1);
        z=cumsum(rk);
        nor=z(end);
        dt=-log(rand(1))/nor;
        frac=rand(1)*nor;
        k=find(z>frac,1);
        X=X+S(:,k);
        t=t+dt;
    end

    tvec=tvec(1:counter-1);
    Xvec=Xvec(:,1:counter-1);
    
    subplot(2,2,i)
    hold on;
    plot(tvec,Xvec(1,:),'DisplayName','$X_1$')
    plot(tvec,Xvec(2,:),'DisplayName','$X_2$')
    plot(tvec,Xvec(3,:),'DisplayName','$X_3$')

    title(sprintf( ' $d_r$ = %.3f  $N$ = %i   ', dr,N));
    ylabel('Population')
    xlabel('Time')
    legend
    hold off;
end

print('1','-dpdf')
