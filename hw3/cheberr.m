% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = cheberr(a,alph,N,v)
    val=sqrt(v/N/(1-alph));
end
