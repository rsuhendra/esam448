% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = clterr(a,alph,N,v)
    val=sqrt(2*v*(erfinv(alph)^2)/N);
end
