clear all, close all;
%%3 SPECTRAL ESTIMATION AND MODELLING

x_1 = randn(128,1);
x_2 = randn(256,1);
x_3 = randn(512,1);  

figure(1)
subplot(3,1,1)
P_1 = pgm(x_1);                 %Test pgm function
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Periodogram of a WGN with 128-sample realisation')
subplot(3,1,2)
P_2 = pgm(x_2);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Periodogram of a WGN with 256-sample realisation')
subplot(3,1,3)
P_3 = pgm(x_3);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Periodogram of a WGN with 512-sample realisation')

%%3.1.1
b = 0.2*[1 1 1 1 1];
P_filt_1 = filter(b,1,P_1);
P_filt_2 = filter(b,1,P_2);
P_filt_3 = filter(b,1,P_3);

figure(2)                   %Filtered periodogram
subplot(3,1,1) 
N=length(P_1);
f=1:1:N;
plot(f/N, P_filt_1);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Filtered periodogram of a WGN with 128-sample realisation')
subplot(3,1,2)
N=length(P_2);
f=1:1:N;
plot(f/N, P_filt_2);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Filtered periodogram of a WGN with 256-sample realisation')
subplot(3,1,3)
N=length(P_3);
f=1:1:N;
plot(f/N, P_filt_3);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Filtered periodogram of a WGN with 512-sample realisation')

%3.1.2
clear x
x = randn(1024,1);
x1 = x(1:128); x2 = x(129:256); x3 = x(257:384); x4 = x(385:512); x5 = x(513:640); x6 = x(641:768); x7 = x(769:896); x8 = x(897:1024);

figure(3)               %128 sample segment periodogram
subplot(2,4,1)
P1 = pgm(x1); 
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,2)
P2 = pgm(x2);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,3)
P3 = pgm(x3);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,4)
P4 = pgm(x4); 
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,5)
P5 = pgm(x5);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,6)
P6 = pgm(x6); 
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,7)
P7 = pgm(x7);
ylabel('Estimated PSD')
xlabel('Normalized frequency')
subplot(2,4,8)
P8 = pgm(x8); 
ylabel('Estimated PSD')
xlabel('Normalized frequency')

Average = (P1+P2+P3+P4+P5+P6+P7+P8)./8;
f=0:1:length(Average)-1;
figure(4)                               %Averaged periodogram
subplot(2,1,1)
plot(f/(length(Average)), Average)
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Averaged periodogram of a 1024-sample sequence of WGN')

Filt_avg = filter(b,1,Average);
subplot(2,1,2)
plot(f/(length(Filt_avg)), Filt_avg)
ylabel('Estimated PSD')
xlabel('Normalized frequency')
title('Filtered average periodogram of a 1024 sample-sequence of WGN')

P1_mean = mean(P1); P1_std = std(P1);
Average_mean = mean(Average); Average_std = std(Average);
Filt_mean = mean(Filt_avg); Filt_std = std(Filt_avg);

%%3.2
clear x
x = randn(1064,1);
b=1;
a=[1 0.9];
y = filter(b, a, x);
y = y(39:length(y));  %intro

figure(5)
time=1:1026;
plot(time, x(39:length(x)))
hold on
plot(time, y)
hold off

[h,w]=freqz([1],[1 0.9],512);           %Theoretical PSD
figure(10)
plot(w/(2*pi),abs(h).^2);
xlabel('w/2\pi ');
ylabel('Squared magnitude response')
hold on
U = pgm(y);


[R lag] = xcorr(y, 'unbiased');
a1 = (-R(find(lag==1)))/(R(find(lag==0)));
var = (R(find(lag==0)))+(a1*R(find(lag==1))); 
[h1,w1] = freqz([var],[1 a1],512);              %Model-based PSD
hold on
plot(w1/(2*pi),abs(h1).^2);     
hold off
legend('Theoretical PSD','Estimated PSD','Model based PSD')


%3.2.5
clear lag
load sunspot.dat
clear y_sun
y_sun = sunspot(:,2);

figure(6)
P = pgm(y_sun);             %sunspot periodogram
f=0:1:length(P)-1;
semilogy(f/length(P),P)     %semilogy used in order to better visualize
[Rx_sun,lag_sun] = xcorr(y_sun,'unbiased');
Rsun = Rx_sun(lag_sun>=0);

%Model based:
p = [1 2 5 10];
 for k=1:length(p)
[a,e] = aryule(y_sun,p(k));
var = a*Rsun(1:p(k)+1);
[h1_orig,w1_orig] = freqz([var],[a],512);
hold on
semilogy(w1_orig/(2*pi),(abs(h1_orig).^2)./length(y_sun));
 end
xlabel('Normalized frequency')
ylabel('PSD')
legend('Periodogram', 'Model based order 1','Model based order 2','Model based order 5','Model based order 10')

%Same analysis for Mean centered sunspot data:
y_mean = sunspot(:,2)-mean(sunspot(:,2));

figure(7)
P = pgm(y_mean);     %Periodogram
f=0:1:length(P)-1;
semilogy(f/length(P),P)

%Model based:
clear a, clear var
p = [1 2 5 10];
 for k=1:length(p)
[a,e] = aryule(y_sun,p(k));
var = a*Rsun(1:p(k)+1);
[h1_orig,w1_orig] = freqz([var],[a],512);
hold on
semilogy(w1_orig/(2*pi),(abs(h1_orig).^2)./length(y_mean));
 end
xlabel('Normalized frequency')
ylabel('PSD')
legend('Periodogram', 'Model based order 1','Model based order 2','Model based order 5','Model based order 10')