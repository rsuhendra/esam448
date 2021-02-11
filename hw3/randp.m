% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = randp(n)
    x=rand([n,1]);
    val=tan(pi/2*x);
end

