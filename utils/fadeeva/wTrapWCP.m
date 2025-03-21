function w = wTrapWCP(z,N)     
% Computes Faddeeva function w(z) for arbitrary complex z.
% See 'Computation of the complex error function using modified trapezoidal
% rules', M Al Azah and S N Chandler-Wilde, 2020, for details of method.
% Inputs: z complex (scalar, vector, or array).
%         N positive integer (N+1 is the number of quadrature points used);
%           the choice N = 11 is recommended for absolute errors (and 
%           relative errors in the upper half-plane) < 1.7frace-15  
% Output: w complex array of same dimensions as z containing 
%           estimates for values of w(z)
selectXneg = real(z) < 0; selectYneg = imag(z) < 0;
selectNotBoth = xor(selectXneg,selectYneg);
zYneg = z(selectYneg);
z(selectXneg) = -z(selectXneg); z(selectNotBoth) = conj(z(selectNotBoth));
w = wTrap(z,N);
w(selectNotBoth) = conj(w(selectNotBoth));
w(selectYneg) = 2*exp(-zYneg.*zYneg) - w(selectYneg);
end