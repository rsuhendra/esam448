function prob = p(x,alph,beta)
    prob = (7/(2*(2+7*beta))).* (x.^2 -x.^20+alph*(x-x.^3)+beta);
end
    