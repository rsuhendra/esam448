% Problem 5 LCG
% (a)
c=1; M=2048; a1=1231; a2=371;
x1 = zeros(1,2); x2 = zeros(1,2);
x1(1)=randi(M); x2(1)=randi(M); t= zeros(1,2);
for n=1:50
    t(n)=n
end
for n = 1:51
    x1(n+1) = mod(a1 * x1(n) +c,M);
    x2(n+1) = mod(a1 * x2(n) +c,M);
end

subplot(2,2,1);
line1=scatter(t,x1(2:51));
ylim([0 M]);
title('xn vs n, a=1231')
subplot(2,2,2); 
line2=scatter(t,x2(2:51));
ylim([0 M]);
title('xn vs n, a=371')
subplot(2,2,3); 
line3=scatter(x1(2:51)/M,x1(3:52)/M);
title('xn+1 vs xn, a=1231 normalized')
subplot(2,2,4); 
line4=scatter(x2(2:51)/M,x2(3:52)/M);
title('xn+1 vs xn, a=371 normalized')

% (b)
N=100; Nt=10000;
co= zeros(1,Nt);
for i= 1:Nt
    x = zeros(1,2);
    x(1)=randi(M);
    for j = 1:N
        x(j+1) = mod(a1 * x(j) +c,M);
    end
    x=x(2:end);
    A=cov(x(1:2:end),x(2:2:end));
    co(i)=A(1,2);
end

subplot(1,1,1);
histogram(co,20);
line([mean(co), mean(co)], ylim, 'LineWidth', 2, 'Color', 'r');
legend(sprintf('stde=%f',std(co)/sqrt(length(co))),sprintf('mean=%f',mean(co)));
title(sprintf('N=%.0f, Nt=%.0f',N,Nt));
