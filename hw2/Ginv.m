% The inverse CDF of the normalized distribution g
function val = Ginv(x,a,b)
    k=(2*(x>0.5)-1);
    val=(k/b).*(-a+sqrt(a^2-k*(2*a*b+b^2)+(k*(4*a*b+2*b^2)).*x));
end

