% CDF of p
function val = bigP(x,alph,beta)
    val = (7/(2*(2+7*beta))).* (x.^3/3 -x.^21/21+alph*(x.^2/2- x.^4/4)+beta*x +1/3-1/21-alph*(1/4)+beta);
end
    