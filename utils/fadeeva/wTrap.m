function w = wTrap(z,N)     
% Computes Faddeeva function w(z) for z = x+iy with x,y >=0.
% See 'Computation of the complex error function using modified trapezoidal
% rules', M Al Azah and S N Chandler-Wilde, 2021, for details of method.
% Inputs: z complex (scalar, vector, or array).
%         N positive integer (N+1 is the number of quadrature points used);
%           the choice N = 11 is recommended for absolute and relative 
%           errors < 2e-15  
% Output: w complex array of same dimensions as z containing 
%           estimates for values of w(z)
h = sqrt(pi/(N+1)); 
H = pi/h;  
rz = real(z); rzh = rz/h; iz = imag(z); 
buff = abs(rzh-floor(rzh)-0.5);
select1 = iz >= max(rz,H); 
select2 = (iz < rz) & (buff <= 0.25);
select3 = ~(select1|select2);
w = zeros(size(z)); 
w(select1) = wM(z(select1),N);
w(select2) = wMT(z(select2),N);
w(select3) = wMM(z(select3),N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ih = wM(z,N) 
% Midpoint rule approximation
z2 = z.*z; az = (2i/H)*z;
t = h*((N:-1:1)+0.5); t2 = t.^2; et2 = exp(-t2); h0 = 0.5*h;
Sum = exp(-h0^2)./(z2-h0^2); 
for n = 1 : N
    Sum = Sum + et2(n)./(z2-t2(n));
end
Ih = az.*Sum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ihstar = wMM(z,N)
% Modified midpoint rule approximation
z2 = z.*z; az = (2i/H)*z;
t = h*((N:-1:1)+0.5); t2 = t.^2; et2 = exp(-t2); h0 = 0.5*h;
Sum = exp(-h0.^2)./(z2-h0.^2); 
for n = 1 : N
    Sum = Sum + et2(n)./(z2-t2(n));
end
Ch = 2./(exp(z2).*(1+exp((-2i*H)*z))); % This is the modification
Ihstar = az.*Sum + Ch;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ihstar = wMT(z,N)   
% Modified trapezium rule approximation
z2 = z.*z; az = (2i/H)*z;
tau = h*(N:-1:1); tau2 = tau.^2; et2 = exp(-tau2);  
Sum = et2(1)./(z2-tau2(1));
for n = 2 : N
    Sum = Sum + et2(n)./(z2-tau2(n));
end
Ch2 = 2./(exp(z2).*(1-exp((-2i*H)*z))); % This is the modification
Ihstar = (1i/H)./z + az.*Sum + Ch2; 
end
end
