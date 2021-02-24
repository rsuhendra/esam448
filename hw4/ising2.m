H=0; J=1; L=20; Nt=500; 
Tc=3:-0.025:1.5;
Lat=sign(rand(L)-0.5);
posit=1:L;
up_shift=circshift(posit,1);
down_shift=circshift(posit,-1);

% Initializing total energy U
U=-H*sum(Lat,'all');
s=0;
for i=1:L
    for j=1:L
        s=s+Lat(i,j)*(Lat(up_shift(i),j)+Lat(down_shift(i),j)+Lat(i,up_shift(j))+Lat(i,down_shift(j)));
    end
end
U=U-(J/2)*s;

% Variables to keep track of
Tcorr=zeros(size(Tc)); Tcorr(1)=20;
SpinFrac=zeros(size(Tc));
SpinEnergy=zeros(size(Tc));
SpecificHeat=zeros(size(Tc));
Magnetization=zeros(size(Tc));
Suscept=zeros(size(Tc));

% Processing
for l=1:length(Tc)
    T=Tc(l);
    tcorr=Tcorr(l);
    Uvec=zeros([Nt*tcorr,1]);
    nrg=zeros([1,Nt]);
    mag=zeros([1,Nt]);
    totspin=0;
    counter=0;
    for z=1:tcorr % Equilibriation step
        [row,col]=ind2sub([L,L],randperm(L^2));
        r=rand(1,L^2);
        for i=1:L^2
            delU= 2*H*Lat(row(i),col(i)) + 2*J*Lat(row(i),col(i))*(Lat(up_shift(row(i)),col(i))+ Lat(down_shift(row(i)),col(i))+ Lat(row(i),up_shift(col(i))) + Lat(row(i),down_shift(col(i))));
            prob=min(1,exp(-delU/T));
            if r(i)<=prob
                Lat(row(i),col(i))=-Lat(row(i),col(i));
                U=U+delU;
            end
        end
    end
    for n=1:Nt
        for z=1:tcorr
            [row,col]=ind2sub([L,L],randperm(L^2));
            r=rand(1,L^2);
            for i=1:L^2
                delU= 2*H*Lat(row(i),col(i)) + 2*J*Lat(row(i),col(i))*(Lat(up_shift(row(i)),col(i))+ Lat(down_shift(row(i)),col(i))+ Lat(row(i),up_shift(col(i))) + Lat(row(i),down_shift(col(i))));
                prob=min(1,exp(-delU/T));
                if r(i)<=prob
                    Lat(row(i),col(i))=-Lat(row(i),col(i));
                    U=U+delU;
                    totspin=totspin+1;
                end
            end
            counter=counter+1;
            Uvec(counter)=U;
        end
        nrg(n)=U;
        mag(n)=abs(sum(Lat,'all'));
    end
    [c,lags]=xcov(Uvec,3000);
    [m,k]=max(c);
    for i=k:length(c)
        if c(i)<(m/10)
            Tcorr(l+1)=i-k;
            break
        end
    end
    SpinFrac(l)=totspin/(Nt*tcorr*L^2);
    SpinEnergy(l)=mean(nrg)/(L^2);
    Magnetization(l)=mean(mag)/(L^2);
    SpecificHeat(l)=(mean(nrg.^2)-mean(nrg)^2)/(T^2);
    Suscept(l)=(mean(mag.^2)-mean(mag)^2)/T;
end
Tcorr(1)=[];

subplot(3,2,1)
x = linspace(0,10);
plot(Tc,SpinEnergy,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Energy per spin')

subplot(3,2,2)
plot(Tc,SpecificHeat,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Specific Heat per spin')

subplot(3,2,3)
plot(Tc,Magnetization,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Magnetization per spin')

subplot(3,2,4)
plot(Tc,Suscept,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Susceptibility per spin')

subplot(3,2,5)
plot(Tc,SpinFrac,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Fraction of accepted spin flips')

subplot(3,2,6)
plot(Tc,Tcorr,'-o')
xlabel('$\tilde{T}$','interpreter','latex')
title('Correlation time')

sgtitle('$L=20, \tilde{T} \in [1.5,3]$','interpreter','latex')
print('FillPageFigure','-dpdf','-fillpage')