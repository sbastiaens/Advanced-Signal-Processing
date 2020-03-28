%3.3 The Least Square estimation of AR coefficients

clear all, close all
%3.3.3
load sunspot.dat
Data = (sunspot(:,2)-mean(sunspot(:,2)))/std(sunspot(:,2));

[rxx, lag] = xcorr(Data,'unbiased');
M = length(Data)-1;
H = zeros(M,10);
x=zeros(M,1);

 for i=1:M
    x(i,1) = rxx(find(lag==i));
end   

for p=1:10
    for i=1:M
    H(i,p) = rxx(find(lag==(i-p)));
    end
end
 a_tot_stand=0;
for p=1:10
a = inv(((H(:,1:p).')*H(:,1:p)))*(H(:,1:p).')*x;
a_tot_stand=[a_tot_stand, a'];
end

