%%%OPTIMAL FILTERING-FIXED AND ADAPTIVE

clear all, close all

%%4.1
%4.1.1
x = randn(1,1000);
a = 1;
b = [1 2 3 2 1];
y = filter(b,a,x)/std(filter(b,a,x)); %to noramlize filter divide by the std
noise = randn(1,1000)*0.1;
z = y + noise;
SNR = 10*log10(var(z)/var(noise));

[rxx lag]= xcorr(x,4,'unbiased');
for i=1:5
    for j=1:5
        Rxx(i,j) = rxx(find(lag==((j-1)*(-1)+(i-1))));  %ACF
    end
end

clear lag
[pzx lag] = xcorr(x,z,4,'unbiased');
for i=1:5
        Pzx(i,1) = pzx(find(lag==((i-1)*(-1))));
end
wopt = inv(Rxx)*Pzx*std(filter(b,a,x)); %need to scale it up by the same amount that y was scaled down


%4.1.2
N(1,:) =  randn(1,1000)*sqrt(0.1);
N(2,:) =  randn(1,1000)*sqrt(0.5);
N(3,:) =  randn(1,1000);
N(4,:) =  randn(1,1000)*sqrt(2);
N(5,:) =  randn(1,1000)*sqrt(5);
N(6,:) =  randn(1,1000)*sqrt(10);

for i=1:6
Z_N(i,:) = y(1,:) + N(i,:);
end
for i=1:6
    SNR_N(i,:) = 10*log10(var(Z_N(i,:))/var(N(i,:))); 
end

for i=1:6
[pzx_N(i,:) lag_N(i,:)] = xcorr(x,Z_N(i,:),4,'unbiased');
    for j=1:5
        Pzx_N(j,i) = pzx_N(i,find(lag==((j-1)*(-1))));
    end
end
wopt_N = inv(Rxx)*Pzx_N*std(filter(b,a,x));

%For Nw greater than 4
[rxx_6 lag]= xcorr(x,6,'unbiased');
for i=1:7
    for j=1:7
        Rxx_6(i,j) = rxx_6(find(lag==((j-1)*(-1)+(i-1))));
    end
end

clear lag
[pzx_6 lag] = xcorr(x,z,6,'unbiased');
for i=1:7
        Pzx_6(i,1) = pzx_6(find(lag==((i-1)*(-1))));
end


[rxx_8 lag]= xcorr(x,8,'unbiased');
for i=1:9
    for j=1:9
        Rxx_8(i,j) = rxx_8(find(lag==((j-1)*(-1)+(i-1))));
    end
end

clear lag
[pzx_8 lag] = xcorr(x,z,8,'unbiased');
for i=1:9
        Pzx_8(i,1) = pzx_8(find(lag==((i-1)*(-1))));
end


[rxx_10 lag]= xcorr(x,10,'unbiased');
for i=1:11
    for j=1:11
        Rxx_10(i,j) = rxx_10(find(lag==((j-1)*(-1)+(i-1))));
    end
end

clear lag
[pzx_10 lag] = xcorr(x,z,10,'unbiased');
for i=1:11
        Pzx_10(i,1) = pzx_10(find(lag==((i-1)*(-1))));
end
wopt_6 = inv(Rxx_6)*Pzx_6*std(filter(b,a,x));
wopt_8 = inv(Rxx_8)*Pzx_8*std(filter(b,a,x));
wopt_10 = inv(Rxx_10)*Pzx_10*std(filter(b,a,x)); %over 5 values are close to 0

%%4.2 LMS algorithm
%for different adaptation gain
[ans,e_lms,w_lms]=lms(x',z',0.001,5);
figure(1)
subplot(3,1,1)
plot(z)
hold on 
plot(ans)
xlabel('Sample number')
ylabel('Signal Amplitude')
title('Signal output z[n] and LMS estimate with \mu = 0.001')
subplot(3,1,2)
[ans,e_lms,w_lms]=lms(x',z',0.01,5);
plot(z)
hold on 
plot(ans)
xlabel('Sample number')
ylabel('Signal Amplitude')
title('Signal output z[n] and LMS estimate with \mu = 0.01')
subplot(3,1,3)
[ans,e_lms,w_lms]=lms(x',z',0.1,5);
plot(z)
hold on 
plot(ans)
xlabel('Sample number')
ylabel('Signal Amplitude')
title('Signal output z[n] and LMS estimate with \mu = 0.1')


%Evolution of coefficients
figure(2)
[ans,e_lms,w_lms]=lms(x',z',0.01,5);
plot(w_lms(1,:)*std(filter(b,a,x))) %scale it up due to y being scaled down
hold on
plot(w_lms(2,:)*std(filter(b,a,x)))
hold on
plot(w_lms(3,:)*std(filter(b,a,x)))
hold on
plot(w_lms(4,:)*std(filter(b,a,x)))
hold on
plot(w_lms(5,:)*std(filter(b,a,x)))
hold off
xlabel('Sample number')
ylabel('Coefficient value')
title('Plot of the 5 coefficients over time estimated with LMS algorithm')
legend('Coefficient 1','Coefficient 2','Coefficient 3','Coefficient 4','Coefficient 5')
figure                      %depending on mu, the solution will reach wiener solution faster or slower.
plot((e_lms).^2)
xlabel('Sample number')
ylabel('Squared estimate error')
title('Plot of the estimated error over time')

[ans,e_lms,w_lms]=lms(x',z',0.002,5);

figure(3)
plot(w_lms(1,:)*std(filter(b,a,x))) 
hold on
plot(w_lms(2,:)*std(filter(b,a,x)))
hold on
plot(w_lms(3,:)*std(filter(b,a,x)))
hold on
plot(w_lms(4,:)*std(filter(b,a,x)))
hold on
plot(w_lms(5,:)*std(filter(b,a,x)))
hold off
xlabel('Sample number')
ylabel('Coefficient value')
legend('Coefficient 1','Coefficient 2','Coefficient 3','Coefficient 4','Coefficient 5')
title('Coefficients of LMS estimate with \mu = 0.002')


[ans,e_lms,w_lms]=lms(x',z',0.1,5);

figure(4)
plot(w_lms(1,:)*std(filter(b,a,x))) %scale it up due to y being scaled down
hold on
plot(w_lms(2,:)*std(filter(b,a,x)))
hold on
plot(w_lms(3,:)*std(filter(b,a,x)))
hold on
plot(w_lms(4,:)*std(filter(b,a,x)))
hold on
plot(w_lms(5,:)*std(filter(b,a,x)))
hold off
xlabel('Sample number')
ylabel('Coefficient value')
legend('Coefficient 1','Coefficient 2','Coefficient 3','Coefficient 4','Coefficient 5')
title('Coefficients of LMS estimate with \mu = 0.1')

[ans,e_lms,w_lms]=lms(x',z',0.2,5);

figure(5)
plot(w_lms(1,:)*std(filter(b,a,x))) %scale it up due to y being scaled down
hold on
plot(w_lms(2,:)*std(filter(b,a,x)))
hold on
plot(w_lms(3,:)*std(filter(b,a,x)))
hold on
plot(w_lms(4,:)*std(filter(b,a,x)))
hold on
plot(w_lms(5,:)*std(filter(b,a,x)))
hold off
xlabel('Sample number')
ylabel('Coefficient value')
legend('Coefficient 1','Coefficient 2','Coefficient 3','Coefficient 4','Coefficient 5')
title('Coefficients of LMS estimate with \mu = 0.2')

%4.3 Gear shifting
[ans,e_lms,w_lms]=lms(x',z',0.1,5); %for gear shifting in the lms function
%mu is varied.

figure(6)
subplot(1,2,1)
plot(w_lms(1,:)*std(filter(b,a,x))) %scale it up due to y being scaled down
hold on
plot(w_lms(2,:)*std(filter(b,a,x)))
hold on
plot(w_lms(3,:)*std(filter(b,a,x)))
hold on
plot(w_lms(4,:)*std(filter(b,a,x)))
hold on
plot(w_lms(5,:)*std(filter(b,a,x)))
hold off
xlabel('Sample number')
ylabel('Coefficient value')
legend('Coefficient 1','Coefficient 2','Coefficient 3','Coefficient 4','Coefficient 5')
title('Coefficients of LMS estimate with varying \mu starting at 0.1 and decreasing by 0.0001 every step')

subplot(1,2,2)
plot((e_lms).^2)
xlabel('Sample number')
ylabel('Squared error')
title('Squared error of the LMS estimate as a function of time')
