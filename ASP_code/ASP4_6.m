%4.6 Dealing with computational complexity: sign algorithms

clear all, close all

w_noise = randn(1,10000);
a = [1 0.9 0.2];
x = filter(1,a,w_noise);

[w, w1, w2, w3] = signalg(x',0.001,3);
figure(1)
subplot(2,2,2)
plot(w1(1,:))
 hold on 
 plot(w1(2,:))
 xlabel('Sample number')
ylabel('Coefficient value')
title('Evolution of a1 and a2 with signed-error LMS algorithm')
subplot(2,2,3)
plot(w2(1,:))
 hold on 
 plot(w2(2,:))
 xlabel('Sample number')
ylabel('Coefficient value')
title('Evolution of a1 and a2 with signed-regressor LMS algorithm')
 subplot(2,2,4)
 plot(w3(1,:))
 hold on 
 plot(w3(2,:))
 xlabel('Sample number')
ylabel('Coefficient value')
title('Evolution of a1 and a2 with sign-sign LMS algorithm')
  subplot(2,2,1)
 plot(w(1,:))
 hold on 
 plot(w(2,:))
xlabel('Sample number')
ylabel('Coefficient value')
title('Evolution of a1 and a2 with basic LMS algorithm')

%Speech comparison
[t,Fs] = audioread('t.m4a');
t_dat = t(52000:53001,1);
[y, y1, y2, y3] = signalgspeech(t_dat,0.1,10);
figure(2)
subplot(2,2,1)
plot(1:length(t_dat),t_dat)
hold on
plot(1:length(y),y)
xlabel('Sample number')
ylabel('Amplitude')
title('Basic LMS algorithm')
legend('Original data', 'Predicted data')
subplot(2,2,2)
plot(1:length(t_dat),t_dat)
hold on
plot(1:length(y),y1)
xlabel('Sample number')
ylabel('Amplitude')
title('Signed-error LMS algorithm')
legend('Original data', 'Predicted data')
subplot(2,2,3)
plot(1:length(t_dat),t_dat)
hold on
plot(1:length(y2),y2)
xlabel('Sample number')
ylabel('Amplitude')
title('Signed-regressor LMS algorithm')
legend('Original data', 'Predicted data')
subplot(2,2,4)
plot(1:length(t_dat),t_dat)
hold on
plot(1:length(y3),y3)
xlabel('Sample number')
ylabel('Amplitude')
title('Sign-sign LMS algorithm')
legend('Original data', 'Predicted data')