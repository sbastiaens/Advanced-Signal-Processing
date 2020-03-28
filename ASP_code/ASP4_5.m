%4.5 Speech recognition

% In comments in enables the importation of the audio file for analysis

%  [e,Fs] = audioread('e.m4a');
%  [a,Fs] = audioread('a.m4a');
%  [s,Fs] = audioread('s.m4a');
%  [t,Fs] = audioread('t.m4a');
 %[xrec,Fs] = audioread('x.m4a');

a_dat = a(79000:80001,1); %To have N=1000, Fs is by default 44100
e_dat = e(56000:57001,1);
s_dat = s(56000:57001,1);
t_dat = t(52000:53001,1);
x_dat = xrec(44000:45001,1);

%Compare the different sounds
figure(1)
subplot(2,3,1)
plot(1:length(a_dat),a_dat(:,1))
[y_a, preda] = speech_rec(a_dat,0.5,10);
hold on
plot(1:length(y_a),y_a)
xlabel('Sample number')
ylabel('Amplitude of signal')
legend('Orginidal data', 'Predicted data')
title('Evolution of data for letter "a"')
xlim([0 1000])
subplot(2,3,2)
plot(1:length(e_dat),e_dat(:,1))
[y_e, prede] = speech_rec(e_dat,0.5,10);
hold on
plot(1:length(y_e),y_e)
xlabel('Sample number')
ylabel('Amplitude of signal')
legend('Orginidal data', 'Predicted data')
title('Evolution of data for letter "e"')
xlim([0 1000])
subplot(2,3,3)
plot(1:length(s_dat),s_dat(:,1))
[y_s, preds] = speech_rec(s_dat,0.5,10);
hold on
plot(1:length(y_s),y_s)
xlabel('Sample number')
ylabel('Amplitude of signal')
legend('Orginidal data', 'Predicted data')
title('Evolution of data for letter "s"')
xlim([0 1000])
subplot(2,3,4)
plot(1:length(t_dat),t_dat(:,1))
[y_t, predt] = speech_rec(t_dat,0.5,10);
hold on
plot(1:length(y_t),y_t)
xlabel('Sample number')
ylabel('Amplitude of signal')
legend('Orginidal data', 'Predicted data')
title('Evolution of data for letter "t"')
xlim([0 1000])
subplot(2,3,[5 6])
plot(1:length(x_dat),x_dat(:,1))
[y_x, predx] = speech_rec(x_dat,0.5,10);
hold on
plot(1:length(y_x),y_x)
xlabel('Sample number')
ylabel('Amplitude of signal')
legend('Orginidal data', 'Predicted data')
title('Evolution of data for letter "x"')
xlim([0 1000])

%Analysis of adaptation gain with letter "t"
figure(2)
subplot(1,3,1)
plot(1:length(t_dat),t_dat(:,1))
hold on
y_t = speech_rec(t_dat,0.5,10);
plot(1:length(y_t),y_t)
xlabel('Sample number')
ylabel('Amplitude of signal')
title('\mu = 0.5')
subplot(1,3,2)
plot(1:length(t_dat),t_dat(:,1))
hold on
y_t = speech_rec(t_dat,5,10);
plot(1:length(y_t),y_t)
xlabel('Sample number')
ylabel('Amplitude of signal')
title('\mu = 5')
subplot(1,3,3)
plot(1:length(t_dat),t_dat(:,1))
hold on
y_t = speech_rec(t_dat,50,10);
plot(1:length(y_t),y_t)
xlabel('Sample number')
ylabel('Amplitude of signal')
title('\mu = 50')