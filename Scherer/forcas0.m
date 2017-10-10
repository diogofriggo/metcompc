fx = zeros(1,N);
fy = zeros(1,N);
fpx = zeros(N,N); fpy = fpx; %forca_x que i faz em j
for k=1:N-1 %1..9
    p = (k+1):N; %2..10, 3..10, ... 9,10
    delx = x(k)-x(p)'; %[x(1) x(1) ... x(1)] - [x(2) x(3)...x(9)] then [x(2) x(2) ... x(2)] - [x(3) x(4)...x(9)] ... [x(9)] - [x(10)]
    dely = y(k)-y(p)';
    R2 = delx.^2 + dely.^2; %[(x(1)-x(2))^2 (x(1)-x(3))^2 ... (x(1)-x(9))^2] up to [(x(8)-x(9))^2 (x(8)-x(10))^2] finally [(x(9)-x(10))^2]
    %matlab is column-major
    list = find(R2<l2); %find indices
    %!isempty(list)
    if(length(list) > 0)
        list = k+list;
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