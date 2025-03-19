% Test of accuracy and correct implementation of wTrapWCP.m by 
% comparing against truncated power series for w(z) given in (2.16) of
% Poppe and Wijers, ACM TOMS (1990), taking N = 12
points = (-200:200)/400; % Uniform grid on [-1/2,1/2] including 0
[x,y] = meshgrid(points); 
% Now set z to be matrix of complex numbers covering {z:x,y\in [-1/2,1/2]}, 
% including points on the positive and negative real and imaginary axes
z = complex(x,y);
wN = wTrapWCP(z,11); % Computing w_N(z) with N = 11
q = z.^2;
w = (((((((((((q/11975040000+1/918086400).*q+1/76204800).*q+ ...
    1/6894720).*q+1/685440).*q+1/75600).*q+1/9360).*q+1/1320).*q+ ...
    1/216).*q+1/42).*q+1/10).*q+1/3).*q+1;
% Now compute w(z) by Poppe and Wijers (2.16) with N = 12
w = exp(-q).*(1+(2i/sqrt(pi))*z.*w); 
ErrorPS = max(max(abs(wN-w))) % Maximum of |w_N(z)-w(z)|