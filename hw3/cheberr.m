% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = cheberr(a,alph,N)
    val=sqrt(varg(a)/N/(1-alph));
end
