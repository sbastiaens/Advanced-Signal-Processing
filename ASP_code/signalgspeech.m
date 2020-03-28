function [y, y_e, y_r, y_s] = signalgspeech(x,mu,order)

w = zeros(order-1,1);
y = zeros(1, length(x));
w_e = zeros(order-1,1);
y_e = zeros(1, length(x));

w_r = zeros(order-1,1);
y_r = zeros(1, length(x));

w_s = zeros(order-1,1);
y_s = zeros(1, length(x));

n=1;
for i=order:length(x)
    xsum = x(i-1:-1:i-order+1);    %basic LMS
    y(i) = w'*xsum;
    e(i) = x(i)-y(i);
    w = w + mu*e(i)*xsum; 
    w_no(:,n) = w;
    
    y_e(i) = w_e'*xsum;             %signed-error
    ee(i) = x(i)-y_e(i);
    w_e = w_e + mu*sign(ee(i))*xsum; 
    w_err(:,n) = w_e;
    
    y_r(i) = w_r'*xsum;             %signed-regressor
    er(i) = x(i)-y_r(i);
    w_r = w_r + mu*er(i)*sign(xsum); 
    w_reg(:,n) = w_r;
    
     y_s(i) = w_s'*xsum;            %sign-sign
    es(i) = x(i)-y_s(i);
    w_s = w_s + mu*sign(es(i))*sign(xsum); 
    w_sign(:,n) = w_s;
    
    n = n+1;
end
end
