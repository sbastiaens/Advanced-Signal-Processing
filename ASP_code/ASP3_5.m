%3.5 Real world signals: Respiratory sinus arrythmia from RR-Intervals

clear all, close all
load('RRI_1')
load('RRI_2')
load('RRI_3')

%Periodogram of the three trials using pgm function
figure(1)                   
subplot(1,3,1)
P1 =  pgm(RR1_1);
f=0:1:length(RR1_1)-1;
plot(f/length(RR1_1), P1)
ylim([0 2])
ylabel('PSD')
xlabel('Normalized frequency')
subplot(1,3,2)
P = pgm(RR1_2);
f=0:1:length(RR1_2)-1;
plot(f/length(RR1_2), P)
ylim([0 2])
ylabel('PSD')
xlabel('Normalized frequency')
subplot(1,3,3)
P3 = pgm(RR1_3);
f=0:1:length(RR1_3)-1;
plot(f/length(RR1_3), P3)
ylim([0 2])
ylabel('PSD')
xlabel('Normalized frequency')


% Averaged periodogram with 50s window
r1 = RR1_1(1:200); r2=RR1_1(201:400); r3 = RR1_1(401:600); r4=RR1_1(601:650);
P1 = pgm(r1);
P2 = pgm(r2);
P3 = pgm(r3);
P4 = pgm(r4);
Average = (P1+P2+P3)./3;
 f=0:1:length(Average)-1;
figure(2)
subplot(1,3,1)
 plot(f/(length(Average)), (Average))
 ylabel('PSD')
xlabel('Normalized frequency')
title('Averaged periodogram trial 1')
 ylim([0 2])


r1 = RR1_2(1:200); r2=RR1_2(201:400); r3 = RR1_2(401:600); r4=RR1_2(601:800);
P1 = pgm(r1);
P2 = pgm(r2);
P3 = pgm(r3);
P4 = pgm(r4);

Average = (P1+P2+P3+P4)./4;
 f=0:1:length(Average)-1;
subplot(1,3,2)
 plot(f/(length(Average)), (Average))
 ylabel('PSD')
xlabel('Normalized frequency')
title('Averaged periodogram trial 2')
 ylim([0 2])


r1 = RR1_3(1:200); r2=RR1_3(201:400); r3 = RR1_3(401:600); r4=RR1_3(601:800);
P1 = pgm(r1);
P2 = pgm(r2);
P3 = pgm(r3);
P4 = pgm(r4);

Average = (P1+P2+P3+P4)./4;
 f=0:1:length(Average)-1;
subplot(1,3,3)
 plot(f/(length(Average)), (Average))
 ylabel('PSD')
xlabel('Normalized frequency')
title('Averaged periodogram trial 3')
 ylim([0 2])

%Note: for a better view a log scale can be chose or zooming in.