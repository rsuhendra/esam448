function rk = reaction(X,N,dr,d,g,a0,a1)
    st=circshift([1 2 3], 1)+3;
    rk=zeros([1,15]);
    for i=1:3
        rk(i)=g*(N-X(st(i)));
        rk(3+i)=d*X(i);
        rk(6+i)=(a0/N)*X(i)*(N-X(3+i));
        rk(9+i)=a1*(X(3+i));
        rk(12+i)=dr*X(3+i);
    end
end