points = linspace(0,10,4001);
[x,y] = meshgrid(points); 
% Now set z to be matrix of complex numbers covering {z:x,y\in [0,10]}, 
% including points on the positive real and imaginary axes
z = complex(x,y);
f0 = @()(exp(-z.*z));
f1 = @()(wTrap(z,11));
f2 = @()(cef(z,40));  
f3 = @()(fadsamp(z));
f4 = @()(Faddeyeva_v2(z,13));

NN = 25
T0 = zeros(1,NN); T1 = T0; T2 = T0; T3 = T0; T4 = T0;
for n = 1:NN
    n
    T0(n) = timeit(f0)
    T1(n) = timeit(f1)
    T2(n) = timeit(f2)
    T3(n) = timeit(f3)
    T4(n) = timeit(f4)
end
T0stats = [length(T0),mean(T0),median(T0),std(T0)]
T1stats = [length(T1),mean(T1),median(T1),std(T1)]
T2stats = [length(T2),mean(T2),median(T2),std(T2)]
T3stats = [length(T3),mean(T3),median(T3),std(T3)]
T4stats = [length(T4),mean(T4),median(T4),std(T4)]