help dinamol0; clear;
par=input('d? X, Y, l, N, dt, tmax, T] \n');
X=par(1); Y=par(2); l=par(3); N=par(4); dt=par(5); tmax=par(6); T=par(7);
%************
% <inicia0> *
%************    
l2=l^2;
S = X*Y;        %88 = 11*8
a = sqrt(S/N);  %sqrt(88/10) = 2.9665
Nx = ceil(X/a); %11/2.9665 = 3.708 ~ 4
Ny = ceil(Y/a); %8/2.9665 = 2.697 ~ 3
a = min(X/Nx,Y/Ny); %min(11/4,8/3) = min(2.75,2.67) = 2.67
x = zeros(1,N);
y = zeros(1,N);
for n=1:N
    z = (n-1)/Nx; %{0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25}
    x(n) = a*(Nx*(z-floor(z))+.5); %1.33333, 4., 6.66667, 9.33333, 1.33333, 4., 6.66667, 9.33333, 1.33333, 4.
    y(n) = a*(ceil(n/Nx)-.5); %1.33333, 1.33333, 1.33333, 1.33333, 4., 4., 4., 4., 6.66667, 6.66667
end
vx = sqrt(T)*randn(1,N);
vy = sqrt(T)*randn(1,N);
%************
% </inicia0> *
%************    

while(length(tmax) > 0)
%*************
% <integra0> *
%*************    
    dt2 = dt*dt; mdt2 = .5*dt2; h=dt/2;
    itmax = ceil(tmax/dt);
    %************
    % <forcas0> *
    %************
    fx = zeros(1,N);
    fy = zeros(1,N);
    fpx = zeros(N,N); %forca_x que i faz em j
    fpy = zeros(N,N); 
    for k=1:N-1 %1..9
        p = (k+1):N; %2..10, 3..10, ... 9,10
        delx = x(k)-x(p); %[x(1) x(1) ... x(1)] - [x(2) x(3)...x(9)] then [x(2) x(2) ... x(2)] - [x(3) x(4)...x(9)] ... [x(9)] - [x(10)]
        dely = y(k)-y(p);
        R2 = delx.^2 + dely.^2; %[(x(1)-x(2))^2 (x(1)-x(3))^2 ... (x(1)-x(9))^2] up to [(x(8)-x(9))^2 (x(8)-x(10))^2] finally [(x(9)-x(10))^2]
        %matlab is column-major
        list = find(R2<l2); %find indices
        %!isempty(list)
        if(length(list) > 0)
            delx = x(k)-x(list); %operates on all neighbours simultaneously
            dely = y(k)-y(list); %[y(1) y(1) ... y(1)] - [y(3) y(7) ... y(10)]
            %[(x(1)-x(3))^2 + (y(1)-y(3))^2      (x(1)-x(7))^2 + (y(1)-y(7))^2   ...   (x(1)-x(10))^2 + (y(1)-y(10))^2]
            r2 = delx.^2 + dely.^2; %[(x(1)-x(3))^2 (x(1)-x(7))^2 ... (x(1)-x(10))^2]
            flj = 12*(1./(r2.^7)-1./(r2.^4)); %ForceLennardJones, r2^7 = r^14 and r2^4 = 8
            %column = column
            fpx(list, k) = (flj.*delx); %force that neighbours impinge on k
            fpy(list, k) = (flj.*dely); %force * distance  [fpy(3,1) fpy(7,1) ... fpy(10,1)] = [flj{1,3} flj{1,7} ... flj{1,10}]*[x(1)-x(3) x(1)-x(7) ... x(1)-x(10)]
            fpx(k, list) = -fpx(list,k)'; %force that k impinges on neighbours
            fpy(k, list) = -fpy(list,k)'; %transpose column to row
        end
    end
    fx = sum(fpx);
    fy = sum(fpy);
    %*************
    % </forcas0> *
    %*************
    for it=1:itmax
        fxi=fx;fyi=fy;
        x = x + vx * dt + mdt2 * fxi;
        y = y + vy * dt + mdt2 * fyi;
        vx = vx + h * (fx + fxi);
        vy = vy + h * (fy + fyi);

    %*************
    % <reflete0> *
    %*************
        list1 = find(x <= 0);
        if(length(list1) > 0)
            vx(list1) = abs(vx(list1));
            x(list1) = 0;
        end

        list2 = find(x >= X);
        if(length(list2) > 0)
            vx(list2) = -abs(vx(list2));
            x(list2) = X;
        end

        list3 = find(y <= 0);
        if(length(list3) > 0)
            vy(list3) = abs(vy(list3));
            y(list3) = 0;
        end

        list4 = find(y >= Y);
        if(length(list4) > 0)
            vy(list4) = -abs(vy(list4));
            y(list4) = Y;
        end
    %**************
    % </reflete0> *
    %**************            

        plot(x, y, 'o')
        pause(0.01)%
    end
%**************
% </integra0> *
%**************        

%         Ec = (sum(vx.^2)+sum(vy))
%         %potencial0;
%         Energia = V+Ec;
%         EcEpE = [Ec, V, Energia]
%         T = Ec/N
%         tmax = input('novo tmax? se "enter", encerra: \n');
%         if(length(tmax) > 0);
%             renorm = sqrt(Tnova/T);
%             vx = vx*renorm;
%             vx = vy*renorm;
%             T = Tnova
%         end
end