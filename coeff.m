N = 50;
beta = 2;

a = zeros(N+1,1);
for i=1:N+1
    a(i) = 1/(factorial(i)*factorial(i-1));
end

b = zeros(N+1,1);
b(1) = 1;
b(2) = 1;
for i=3:N+1
    b(i) = (b(i-1) - (2*i-3)*a(i-1))/((i-1)*(i-2));
end

dx = .01;
x = 0:dx:1;

y1 = a(1);
for i=1:N
    y1 = y1 + a(i+1)*x.^i;
end
y1p = a(2);
for i=2:N
    y1p = y1p + a(i+1)*i*x.^(i-1);
end
y1pp = 2*a(3);
for i=3:N
    y1pp = y1pp + a(i+1)*i*(i-1)*x.^(i-2);
end

test1 = (x.^2).*y1pp + (2*x).*y1p - (x).*y1;

y2 = y1.*log(x);
for i=0:N
    y2 = y2 + b(i+1)*x.^(i-1);
end

y2p = y1p.*log(x) + y1./x - b(1)./(x.^2);
for i=1:N
    y2p = y2p + b(i+1)*(i-1)*x.^(i-2);
end

y2pp =y1pp.*log(x) + 2*y1p./x - y1./(x.^2) + 2*b(1)./(x.^3);
for i=2:N
    y2pp = y2pp + b(i+1)*(i-2)*(i-1)*x.^(i-3);
end

test2 = (x.^2).*y2pp + (2*x).*y2p - (x).*y2;

f = x.^beta;

f1 = (x.^2).*y2.*f./(y1p.*y2+y1.*y2p);
F1 = zeros(size(x));
F1(1) = 0;
for i=2:length(x)
    F1(i) = dot(f1(2:i),dx*ones(i-1,1));
end

f2 = (x.^2).*y1.*f./(y1p.*y2+y1.*y2p);
F2 = zeros(size(x));
F2(1) = 0;
for i=2:length(x)
    F2(i) = dot(f2(2:i),dx*ones(i-1,1));
end

C1 = (y1p(end)*F1(end) + y2p(end)*F2(end) + y1(end)*f1(end) + y2(end)*f2(end))/y1p(end);

y = y1.*(F1 - C1) +y2.*F2;

figure, clf
subplot(2,2,1)
plot(x,y1)
title('homogeneous solution 1')
subplot(2,2,2)
plot(x,y2)
title('homogeneous solution 2')
subplot(2,2,3)
plot(x,test1)
title('homogeneous solution 1 test')
subplot(2,2,4)
plot(x,test2)
title('homogeneous solution 2 test')
figure, clf
plot(x,y)
title('particular solution')