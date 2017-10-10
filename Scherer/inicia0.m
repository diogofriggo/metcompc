par=input('de X, Y, l, N, dt, tmax, T] \n');
X=par(1); Y=par(2); l=par(3); N=par(4); 
dt=par(5); tmax=par(6); T=par(7); l2 = l^2;
S = X*Y;        %88 = 11*8
a = sqrt(S/N);  %sqrt(88/10) = 2.9665
Nx = ceil(X/a); %11/2.9665 = 3.708 ~ 4
Ny = ceil(Y/a); %8/2.9665 = 2.697 ~ 3
a = min(X/Nx,Y/Ny); %min(11/4,8/3) = min(2.75,2.67) = 2.67
x = zeros(1,N); y = zeros(1,N);
for n=1:N
    z = (n-1)/Nx; %{0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25}
    x(n) = a*(Nx*(z-floor(z))+.5); %1.33333, 4., 6.66667, 9.33333, 1.33333, 4., 6.66667, 9.33333, 1.33333, 4.
    y(n) = a*(ceil(n/Nx)-.5); %1.33333, 1.33333, 1.33333, 1.33333, 4., 4., 4., 4., 6.66667, 6.66667
end
vx = sqrt(T)*randn(1,N);
vy = sqrt(T)*randn(1,N);

