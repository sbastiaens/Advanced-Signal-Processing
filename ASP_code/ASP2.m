clear all, clc, close all
%2 LINEAR STOCHASTIC MODELLING

%%2.1
x = randn(1000,1);
[ACF lag] = xcorr(x, 'unbiased'); %ACFof x
figure(1)
plot(lag, ACF)
axis([-999 999 -0.8 1.2]);
xlabel('Lags')
ylabel('Autocorrealtion')

%4
y=filter(ones(9,1),[1],x);
figure(2)
[xc,lags] = xcorr(y,20,'unbiased');  %ACF of y
stem(lags,xc,'filled')
xlabel('Lags')
ylabel('Autocorrealtion')

y_1=filter(ones(3,1),[1],x); %filter with different order
y_2=filter(ones(5,1),[1],x);
y_3=filter(ones(10,1),[1],x);
y_4=filter(ones(15,1),[1],x);
figure(3)
[xc,lags] = xcorr(y_1,20,'unbiased');  
stem(lags,xc,'filled')
hold on
[xc,lags] = xcorr(y_2,20,'unbiased');
stem(lags,xc,'filled')
hold on
[xc,lags] = xcorr(y_3,20,'unbiased');
stem(lags,xc,'filled')
hold on
[xc,lags] = xcorr(y_4,20,'unbiased');
stem(lags,xc,'filled')
legend('MA = 3','MA = 5','MA = 10','MA = 15')
xlabel('Lags')
ylabel('Autocorrealtion')
hold off

%%2.2
figure(4)
[xcc,lag_cc] = xcorr(x,y,20,'unbiased');
stem(lag_cc,xcc)
ylabel('Cross-correlation of x and y')
xlabel('Lags')
% To compute for other order models to assess the effect
% figure(15)
% [xcc,lag_cc] = xcorr(x,y_1,20,'unbiased');
% stem(lag_cc,xcc)
% ylabel('Cross-correlation of x and y')
% xlabel('Lags')
% figure
% [xcc,lag_cc] = xcorr(x,y_2,20,'unbiased');
% stem(lag_cc,xcc)
% figure
% [xcc,lag_cc] = xcorr(x,y_3,20,'unbiased');
% stem(lag_cc,xcc)
% figure
% [xcc,lag_cc] = xcorr(x,y_4,20,'unbiased');
% stem(lag_cc,xcc)
% ylabel('Cross-correlation of x and y')
% xlabel('Lags')


%%2.3
a1 = -2.5 + 5*rand(100,1);
a2 = -1.5 + 3*rand(100,1);
for i=1:100
 if isstable(1, [1 -a1(i) -a2(i)]); %determines when system is stable 
     figure(5)
     hold on
     plot(a1(i), a2(i), '*')  
 end
end
xlabel('a1')
ylabel('a2')
t = linspace(-2,2); %plot stability triangle
t1(1:100) = -1;
t2 = 0:0.1:2;
t3 = -t2+1;
t4 = -2:0.1:0;
t5 = t4+1;
plot(t,t1,'r')
leg = plot(t,t1,'r')
hold on
plot(t2,t3,'r')
hold on
plot(t4,t5,'r')
legend([leg],{'Stability triangle'})
hold off

% The following part commented is to assess the stability triangle

% a1 = -2.5 + 5*rand(10000,1);
% a2 = -1.5 + 3*rand(10000,1);
% for i=1:10000
%  if isstable(1, [1 -a1(i) -a2(i)]);
%  figure(15)
%      hold on
%     plot(a1(i), a2(i),'r*') %for stability depending on a1 and a2
%      p1 = plot(a1(i), a2(i),'r*') ;
%  else
%   figure(15)
%      hold on
%    plot(a1(i), a2(i),'b*')
%    p2 =  plot(a1(i), a2(i),'b*');
%  end
% end
% xlabel('a1')
% ylabel('a2')
% legend([p1,p2],'\color{red} stable', '\color{blue} unstable')


figure(6)
load sunspot.dat %2nd column is the data
[ACF_sun lag_sun] = xcorr(sunspot(1:5,2), 'unbiased');
subplot(3,1,1)
plot(lag_sun, ACF_sun)
hold on
[ACF_sun lag_sun] = xcorr(sunspot(1:5,2)-mean(sunspot(1:5,2)), 'unbiased');
plot(lag_sun, ACF_sun)
hold off
title('Data length = 5')
ylabel('Autocorrelation')
xlabel('Lag')
[ACF_sun lag_sun] = xcorr(sunspot(1:20,2), 'unbiased');
subplot(3,1,2)
plot(lag_sun, ACF_sun)
hold on
[ACF_sun lag_sun] = xcorr(sunspot(1:20,2)-mean(sunspot(1:20,2)), 'unbiased');
plot(lag_sun, ACF_sun)
hold off
title('Data length = 20')
ylabel('Autocorrelation')
xlabel('Lag')
[ACF_sun lag_sun] = xcorr(sunspot(1:250,2), 'unbiased');
subplot(3,1,3)
plot(lag_sun, ACF_sun)
hold on
[ACF_sun lag_sun] = xcorr(sunspot(1:250,2)-mean(sunspot(1:250,2)), 'unbiased');
plot(lag_sun, ACF_sun)
hold off
title('Data length = 250')
ylabel('Autocorrelation')
xlabel('Lag')
legend('With mean', 'With zero mean')

figure(8)
%Yule Walker
a_tot = 0;
for p=1:10
    a_sun=aryule(sunspot(:,2),p);
    a_tot=[a_tot, a_sun];
end
[~,~,reflect_coeffs] = aryule(sunspot(:,2),10);
stem(1:10, -reflect_coeffs);
axis([0 10.5 -1 1])
title('Partial Autocorrelation function')
ylabel('Partial Autocorrelations')
xlabel('Lag')

%zero mean and unit variance
figure(9)
Sunspot_0_1 = (sunspot(:,2)-mean(sunspot(:,2)))/std(sunspot(:,2)); 
%Zero mean and unit standard deviation
a_tot_stand = 0;
for p=1:10
    a_sun=aryule(Sunspot_0_1,p);
    a_tot_stand=[a_tot_stand, a_sun];
end
[~,~,reflect_coeffs_stand] = aryule(Sunspot_0_1,10);
stem(1:10, -reflect_coeffs_stand);
axis([0 10.5 -1 1])
title('Partial Autocorrelation function')
ylabel('Partial Autocorrelations')
xlabel('Lag')

%MDL AND AIC

N = length(Sunspot_0_1);
%MDL
for p=1:10
 [a_mdl,e] = aryule(Sunspot_0_1,p);
 e_n(p) = e;
MDL(p) =log10(e_n(p))+((p*log10(N))/(N));
AIC(p) =log10(e_n(p))+((2*p)/(N));
AICC(p) = AIC(p)+((2*p*(p+1))/(N-p-1));
end
figure(10)
subplot(3,1,1)
plot(1:10,log10(e_n))
hold on
plot(1:10,MDL)
title('MDL for Standardised Sunspot time series')
xlabel('Model order p')
legend('Cumulative squared error','MDL')
subplot(3,1,2)
plot(1:10,log10(e_n))
hold on
plot(1:10,AIC)
title('AIC for Standardised Sunspot time series')
xlabel('Model order p')
legend('Cumulative squared error','AIC')
subplot(3,1,3)
plot(1:10,log10(e_n))
hold on
plot(1:10,AICC)
title('AICc for Standardised Sunspot time series')
xlabel('Model order p')
legend('Cumulative squared error','AICc')
hold off

%This commented part is to do 5 in order to compare prediction horizons and model order 

% %for AR(1)
% N = length(Sunspot_0_1);
% p = [1, 2, 5, 10];
% for k=1:4
%     for i=p(k)+1:N+1
%     x_pred(k,i+p(k)) = 0.8212*Sunspot_0_1(i-1);
%     end
%     figure
% plot(1:N,Sunspot_0_1)
% hold on
% plot(1:N, x_pred(k,1:288))
% end
% xlabel('Sample number')
% ylabel('Amplitude')
% legend('Actual data','Predicted data')
% a = aryule(Sunspot_0_1,2);
% for k=1:4
%     for i=p(k)+2:N+1
%     x_pred2(k,i+p(k)) = -a(2)*Sunspot_0_1(i-1)-a(3)*Sunspot_0_1(i-2);
%     end
%     figure
% plot(1:N,Sunspot_0_1)
% hold on
% plot(1:N, x_pred2(k,1:288))
% end
% xlabel('Sample number')
% ylabel('Amplitude')
% legend('Actual data','Predicted data')
% 
% a = aryule(Sunspot_0_1,10);
% for k=1:4
%     for i=p(k)+10:N+1
%     x_pred10(k,i+p(k)) = -a(2)*Sunspot_0_1(i-1)-a(3)*Sunspot_0_1(i-2)-a(4)*Sunspot_0_1(i-3)-a(5)*Sunspot_0_1(i-4)-a(6)*Sunspot_0_1(i-5)-a(7)*Sunspot_0_1(i-6)-a(8)*Sunspot_0_1(i-7)-a(9)*Sunspot_0_1(i-8)-a(10)*Sunspot_0_1(i-9)-a(11)*Sunspot_0_1(i-10);
%     end
%     figure
% %plot(1:N,Sunspot_0_1)
% hold on
% plot(1:N, x_pred10(k,1:288))
% end
% xlabel('Sample number')
% ylabel('Amplitude')
% legend('Actual data','Predicted data')


%%2.4
load('NASDAQ.mat')
Date = NASDAQ.Date; 
Data = NASDAQ.Close;
figure(11)
N = length(Data);
%a
parcorr(Data) %goes down after lag 1 so AR(1)
% MDL and AIC
for p=1:10
 [a_mdl,e] = aryule(Data,p);
 e_n(p) = e;
MDL(p) =log10(e_n(p))+((p*log10(N))/(N));
AIC(p) =log10(e_n(p))+((2*p)/(N));
AICC(p) = AIC(p)+((2*p*(p+1))/(N-p-1));
end
figure(12)
subplot(3,1,1)
plot(1:10,log10(e_n))
hold on
plot(1:10,MDL)
title('MDL for AR(p)')
xlabel('Model order')
legend('Cumulative squared error','MDL')
subplot(3,1,2)
plot(1:10,log10(e_n))
hold on
plot(1:10,AIC)
title('AIC for AR(p)')
xlabel('Model order')
legend('Cumulative squared error','AIC')
subplot(3,1,3)
plot(1:10,log10(e_n))
hold on
plot(1:10,AICC)
title('AICc for AR(p)')
xlabel('Model order')
legend('Cumulative squared error','AICc')
hold off

%c)
%i)
var_noise = 1:50:1001;
N = 1:50:1001;
for i=1:21
    for j=1:21
        CRLB_var(i,j) = (2*(var_noise(i)).^(2))/N(j);
    end
end

a1 = -1:0.1:1;
N = 1:50:1001;
for i=1:21
    for j=1:21
        CRLB_a1(i,j) = (1./N(j)).*(1-(a1(i)).^(2));
    end
end
figure(13)
h = heatmap(N,var_noise,CRLB_var,'ColorScaling','log');
h.Title = 'Cramer-Rao lower bound for variance';
h.YLabel = 'True variance';
h.XLabel = 'N';
figure(14)
h = heatmap(N,a1,CRLB_a1,'ColorScaling','log');
h.Title = 'Cramer-Rao lower bound for coefficient a1';
h.YLabel = 'True coefficent a1';
h.XLabel = 'N';



%ii
a_var = aryule(Data,1);
var_Dat = (1/length(Data))*(1+(a_var(2))^2);

%%2.5

%The following comments are to load the ECG trials

% load('RAW.mat')
% Trial_1 = data(find(time==74.43):find(time==237.1));
% Trial_2 = data(find(time==245.2):find(time==489.5));
% Trial_3 = data(find(time==494.3):find(time==736.9));
% [RR1_1, RRIf_1] = ECG_to_RRI(Trial_1,1000);
% [RR1_2, RRIf_2] = ECG_to_RRI(Trial_2,1000);
% [RR1_3, RRIf_3] = ECG_to_RRI(Trial_3,1000);
% save('RRI_1.mat','RR1_1','RRIf_1');
% save('RRI_2.mat','RR1_2','RRIf_2');
% save('RRI_3.mat','RR1_3','RRIf_3');
load('RRI_1')
load('RRI_2')
load('RRI_3')
t=1/RRIf_1:1/RRIf_1:length(RR1_1)/RRIf_1;
figure(15)
plot(t,RR1_1)
load('RRI_3')
t=1/RRIf_2:1/RRIf_2:length(RR1_2)/RRIf_2;
figure(16)
plot(t,RR1_2)
t=1/RRIf_3:1/RRIf_3:length(RR1_3)/RRIf_3;
figure(17)
plot(t,RR1_3)
%RR1_1: data trial 1; RR1_2: data trial 2; RR1_3: data trial 3; 

%heart rate PDE analysis
h_1 = 60./RR1_1;
alpha =1;
for i=1:length(RR1_1)/10
   h_hat_1(i) = (1/10)*sum(alpha*h_1((10*(i-1)+1):10*i)) ;
end
alpha = 0.6;
for i=1:length(RR1_1)/10
   h_hat_2(i) = (1/10)*sum(alpha*h_1((10*(i-1)+1):10*i)); 
end
figure(18)
histogram(h_1,10)
hold on
histogram(h_hat_1)
legend('Original HR', 'Average HR, \alpha = 1')
ylabel('PDE')
xlabel('Heart rate (beats/min)')
hold off
figure(19)
histogram(h_1,10)
hold on
histogram(h_hat_2)
ylabel('PDE')
xlabel('Heart rate (beats/min)')
legend('Original HR', 'Average HR, \alpha = 0.6')
hold off

%ACF of different trials
figure(20)
subplot(3,1,1)
rr1 = detrend(RR1_1);
[ACF_rr1 lag] = xcorr(rr1,'unbiased');
stem(lag, ACF_rr1)
xlim([-400 400])
ylabel('Autocorrelation')
xlabel('Lag')
legend('Trial 1: Normal Breathing')
subplot(3,1,2)
rr2 = detrend(RR1_2);
[ACF_rr2 lag] = xcorr(rr2,'unbiased');
stem(lag, ACF_rr2)
xlim([-400 400])
ylabel('Autocorrelation')
xlabel('Lag')
legend('Trial 2: Fast Breathing')
rr3 = detrend(RR1_3);
[ACF_rr3 lag] = xcorr(rr3,'unbiased');
subplot(3,1,3)
stem(lag, ACF_rr3)
xlim([-400 400])
ylabel('Autocorrelation')
xlabel('Lag')
legend('Trial 3: Deep Breathing')
hold off

%Model order analysis of the three trials
figure(21)
parcorr(rr3)
N=length(rr3);
for p=1:10
 [a_mdl,e] = aryule(rr3,p);
 e_n(p) = e;
MDL(p) =log10(e_n(p))+((p*log10(N))/(N));
AIC(p) =log10(e_n(p))+((2*p)/(N));
AICC(p) = AIC(p)+((2*p*(p+1))/(N-p-1));
end
figure(22)
subplot(3,1,1)
plot(1:10,log10(e_n))
hold on
plot(1:10,MDL)
title('MDL for RRI 3')
xlabel('Model order')
legend('Cumulative squared error','MDL')
subplot(3,1,2)
plot(1:10,log10(e_n))
hold on
plot(1:10,AIC)
title('AIC for RRI 3')
xlabel('Model order')
legend('Cumulative squared error','AIC')
subplot(3,1,3)
plot(1:10,log10(e_n))
hold on
plot(1:10,AICC)
title('AICc for RRI 3')
xlabel('Model order')
legend('Cumulative squared error','AICc')
hold off

figure(23)
parcorr(rr1)

N=length(rr1);
for p=1:10
 [a_mdl,e] = aryule(rr1,p);
 e_n(p) = e;
MDL(p) =log10(e_n(p))+((p*log10(N))/(N));
AIC(p) =log10(e_n(p))+((2*p)/(N));
AICC(p) = AIC(p)+((2*p*(p+1))/(N-p-1));
end
figure(24)
subplot(3,1,1)
plot(1:10,log10(e_n))
hold on
plot(1:10,MDL)
title('MDL for RRI 1')
xlabel('Model order')
legend('Cumulative squared error','MDL')
subplot(3,1,2)
plot(1:10,log10(e_n))
hold on
plot(1:10,AIC)
title('AIC for RRI 1')
xlabel('Model order')
legend('Cumulative squared error','AIC')
subplot(3,1,3)
plot(1:10,log10(e_n))
hold on
plot(1:10,AICC)
title('AICc for RRI 1')
xlabel('Model order')
legend('Cumulative squared error','AICc')
hold off

figure(25)
parcorr(rr2)

N=length(rr2);
for p=1:10
 [a_mdl,e] = aryule(rr2,p);
 e_n(p) = e;
MDL(p) =log10(e_n(p))+((p*log10(N))/(N));
AIC(p) =log10(e_n(p))+((2*p)/(N));
AICC(p) = AIC(p)+((2*p*(p+1))/(N-p-1));
end
figure(26)
subplot(3,1,1)
plot(1:10,log10(e_n))
hold on
plot(1:10,MDL)
title('MDL for RRI 2')
xlabel('Model order')
legend('Cumulative squared error','MDL')
subplot(3,1,2)
plot(1:10,log10(e_n))
hold on
plot(1:10,AIC)
title('AIC for RRI 2')
xlabel('Model order')
legend('Cumulative squared error','AIC')
subplot(3,1,3)
plot(1:10,log10(e_n))
hold on
plot(1:10,AICC)
title('AICc for RRI 2')
xlabel('Model order')
legend('Cumulative squared error','AICc')
hold off