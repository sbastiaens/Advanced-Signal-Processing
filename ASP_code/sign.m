function [w_err, w_reg, w_sign] = sign(x,mu,order)
 
w=zeros(order,N);
y = zeros(1, length(x));
 n=1;
    for i=order:length(x)
       xsum = x(i:-1:i-order+1);
       y(i) = w'*xsum; 
       e(i) = x(i) - y(i);
       we = w + mu*sign(e(i))*xsum; 
       w_err(:,n) = we;
       wr = w + mu*e(i)*sign(xsum); 
       w_reg(:,n) = wr;
       ws = w + mu*sign(e(i))*sign(xsum); 
       w_sign(:,n) = ws;
       n = n+1;
    end
end