function [h,y1] = VDF_polynomial_coefficient(g,N,x)
y = g(x);
h = polyfit(x,y,N);
y1 = zeros(size(y));
for i = 1:N+1
    y1 = h(i)*x.^(N+1-i)+y1;
end
figure;plot(x,y);hold on;plot(x,y1);title('Polynomial transfer function')