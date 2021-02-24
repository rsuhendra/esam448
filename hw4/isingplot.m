H=0; J=1; L=100; Nt=500; 
Tc=[3 2.27 1.5];
Tcorr=[50 1000 50];
Lat=sign(rand(L)-0.5);
posit=1:L;
up_shift=circshift(posit,1);
down_shift=circshift(posit,-1);

% Initializing U
U=-H*sum(Lat,'all');
s=0;
for i=1:L
    for j=1:L
        s=s+Lat(i,j)*(Lat(up_shift(i),j)+Lat(down_shift(i),j)+Lat(i,up_shift(j))+Lat(i,down_shift(j)));
    end
end
U=U-(J/2)*s;

for l=1:length(Tc)
    T=Tc(l);
    tcorr=Tcorr(l);
    for n=1:Nt+1
        for z=1:tcorr
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
    end
    subplot(2,2,4-l)
    pcolor(Lat)
    msg=sprintf('T=%.2f',T);
    title(msg)
end