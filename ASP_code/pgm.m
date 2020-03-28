function [P] = pgm(x)

N = length(x);
f = 0:1:N;
%f = f/N;

for t=1:length(f)
    Abs = 0;
   for j=1:N
     Abs = (x(j)*exp(-1i*2*pi*f(t)*(j/N)))+ Abs;
   end
    P(t) = (1/N)*((abs(Abs))^2);
end
plot(f/N, P)
end