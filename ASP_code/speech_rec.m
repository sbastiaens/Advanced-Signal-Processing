function[y,pred] = speech_rec(x_val,mu,order)
x = x_val;
w_1 = zeros(1,1);
w_2 = zeros(1,1);
y_1 = zeros(1, length(x));
y_2 = zeros(1,length(x));
n=1;
for i=order:length(x)
    x_1 = x(i-1);
    x_2 = x(i-2);
    y_1(i) = w_1'*x_1;
    y_2(i) = w_2'*x_2;
    e(i) = x(i)-(y_1(i)+y_2(i));
    w_1 = w_1 + mu*e(i)*x_1; 
    w_2 = w_2 + mu*e(i)*x_2; 
    w_ad1(:,n) = w_1;
    w_ad2(:,n) = w_2;
    n = n+1;
    y(i) = y_1(i)+y_2(i);
end
pred = 10*log10(var(x)/var(e));
end
