%4.4 Identificqtion of AR processes

clear all, close all;
w_noise = randn(1,1000);
a = [1 0.9 0.2];
x = filter(1,a,w_noise);
w_1 = zeros(1,1);
w_2 = zeros(1,1);
y_1 = zeros(1, length(x));
y_2 = zeros(1,length(x));
order = 3;

%For different adaptation gain
mu = 0.01;
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
end

subplot(2,2,1)
plot(w_ad1)
hold on
plot(w_ad2)
ylabel('Coefficient value')
xlabel('Sample number')
title('Evolution of a1 and a2 over time for \mu = 0.01')

mu = 0.02;
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
end
subplot(2,2,2)
plot(w_ad1)
hold on
plot(w_ad2)
ylabel('Coefficient value')
xlabel('Sample number')
title('Evolution of a1 and a2 over time for \mu = 0.02')

mu = 0.06;
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
end
subplot(2,2,3)
plot(w_ad1)
hold on
plot(w_ad2)
ylabel('Coefficient value')
xlabel('Sample number')
title('Evolution of a1 and a2 over time for \mu = 0.06')

mu = 0.1;
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
end
subplot(2,2,4)
plot(w_ad1)
hold on
plot(w_ad2)
ylabel('Coefficient value')
xlabel('Sample number')
title('Evolution of a1 and a2 over time for \mu = 0.1')