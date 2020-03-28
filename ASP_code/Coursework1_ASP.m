clear all, clc, close all

%1.1
x = rand(1000,1);
figure(1)
plot(1:1000, x)

m = mean(x);
standard_dev = std(x);

for i=1:10
  x = rand(1000,1);
  x_n{i,1} = x(:,1);
end

for i=1:10
    m_n(i) = mean(x_n{i,1});
    std_n(i) = std(x_n{i,1});
end
figure(2)
scatter(1:10, m_n)
hold on
scatter(1:10, std_n)
title('Sample mean and standard deviation for 10 realisations')
xlabel('Realisation number')
ylabel('AU')
legend('sample mean', 'standard deviation')
print -depsc myfig.eps
t=1:0.01:30;

    
for i=1:10
   max(i) = (1/sqrt(12)) - std_n(i) 
end
%hist(x)

%1.1-5
x = randn(1000,1);

m5 = mean(x);
standard_dev5 = std(x);

for i=1:10
  x5 = randn(1000,1);
  x_n5{i,1} = x5(:,1);
end

for i=1:10
    m_n5(i) = mean(x_n5{i,1});
    std_n5(i) = std(x_n5{i,1});
end
figure(3)
scatter(1:10, m_n5)
hold on
scatter(1:10, std_n5)
title('Sample mean and standard deviation using randn')
xlabel('Realisation number')
ylabel('AU')
legend('sample mean', 'standard deviation')
print -depsc myfig2.eps
t=1:0.01:30;

    
for i=1:10
   max(i) = (1/sqrt(12)) - std_n(i) 
end
%hist(x)

