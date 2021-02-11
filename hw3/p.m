% Generates n random numbers distributed 2/pi*1/(1+x^2) x>=0
function val = p(x)
    val=(2/pi).*(1./(1+x.^2));
end