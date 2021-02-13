% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = clterr(a,alph,N)
    val=sqrt(2*varg(a)*(erfinv(alph)^2)/N);
end
