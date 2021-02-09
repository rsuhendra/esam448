function y = reject(alph,beta,n)
    a=0.2; b=1;
    xhat=Ginv(rand(1,n),a,b);
    reject_test=p(xhat,alph,beta)./f(xhat,a,b); 
    uniform=rand(size(xhat)); 
    accept=uniform<=reject_test; 
    y=xhat(accept);
end