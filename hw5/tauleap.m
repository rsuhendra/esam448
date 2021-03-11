N=10000; dr=0.00; d=0.003; g=0.05; a0=0.5; a1=0.01; 
T=15000; tau=0.2;

S=zeros([6,15]); 
for i=1:3
    S(i,i)=1;
    S(i,3+i)=-1;
    S(i,6+i)=-1; S(3+i,6+i)=1;
    S(i,9+i)=1; S(3+i,9+i)=-1;
    S(3+i,12+i)=-1;
end

Xvec=NaN([3,T/tau+2]);
tvec=NaN([1,T/tau+2]);
t=0; X=zeros([6,1]);
counter=1;
while t<T
    tvec(counter)=t;
    Xvec(:,counter)=X(1:3);
    counter=counter+1;
    
    rk=reaction(X,N,dr,d,g,a0,a1);
    poiz=poissrnd(tau*rk);
    
    if any(rk<0)
        disp('stop')
        break
    end
    
    Xmid=X+0.5*tau*S*poiz.';
    rk=reaction(Xmid,N,dr,d,g,a0,a1);
    poiz=poissrnd(tau*rk);
    X=X+S*poiz.';
    
    if any(X<0)
        disp('stop')
        break
    end
    if any(rk<0)
        disp('stop')
        break
    end
    t=t+tau;
end

tvec=tvec(1:counter-1);
Xvec=Xvec(:,1:counter-1);

msX1=Xvec(1,:)-tau*sum(Xvec(1,: ))/T;
X1hat=fft(msX1);
L=length(tvec);
P2 = abs(X1hat/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (2*pi/T)*(0:(L/2));

hold on;
plot(f(1:100),P1(1:100),'DisplayName','$\hat{X_1}$')
legend
title(sprintf( 'Fourier Transform of Mean Subtracted $X_1$'));
ylabel('Power spectrum')
xlabel('$\omega$ (Hz)')
hold off;