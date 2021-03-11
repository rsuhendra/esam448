set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

Ns=[50 100 500]; drs=0:0.0001:0.003; d=0.003; g=0.05; a0=0.5; a1=0.01; 
T=20000;
V=zeros([3,length(drs)]);

S=zeros([6,15]); 
for i=1:3
    S(i,i)=1;
    S(i,3+i)=-1;
    S(i,6+i)=-1; S(3+i,6+i)=1;
    S(i,9+i)=1; S(3+i,9+i)=-1;
    S(3+i,12+i)=-1;
end
for j=1:length(Ns)
    N=Ns(j);
    for i=1:length(drs)
        dr=drs(i);
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

        V(1,i)=var(Xvec(1,:));
        V(2,i)=var(Xvec(2,:));
        V(3,i)=var(Xvec(3,:)); 
    end
    subplot(length(Ns),1,j)
    hold on;
    plot(drs,V(1,:),'-o','DisplayName','$X_1$')
    plot(drs,V(2,:),'-o','DisplayName','$X_2$')
    plot(drs,V(3,:),'-o','DisplayName','$X_3$')
    title(sprintf( '$N$ = %i',N));
    ylabel('$\mathcal{V}(d_r)$')
    xlabel('$d_r$')
    legend
    hold off;

end

print('2','-dpdf')