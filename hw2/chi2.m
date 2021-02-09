function val = chi2(pj,Nj)
    N=sum(Nj);
    Ej=N.*pj;
    vec=(Nj-Ej).^2./Ej;
    val=sum(vec);
end
    