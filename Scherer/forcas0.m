function forcas0(l)
    fx = zeros(1,N);
    fy = zeros(1,N);
    fpx = zeros(N,N); %forca_x que i faz em j
    fpy = zeros(N,N); 
    for k=1:N-1 %1..9
        p = (k+1):N; %2..10, 3..10, ... 9..10
        delx = x(k)-x(p);
        dely = y(k)-y(p);
        R2 = delx.^2 + dely.^2;
        %matlab is column-major
        list = find(R2<l^2); %find indices
        %!isempty(list)
        if(length(list) > 0)
            delx = x(k)-x(list); %operates on all neighbours simultaneously
            dely = x(k)-x(list);
            r2 = delx.^2+dely.^2;
            flj = 12*(1./(r2.^7)-1./(r2.^4)); %fORCElENNARDjONES, r2^7 = r^14 and r2^4 = 8
            %column = column
            fpx(list, k) = (flj.*delx); %force that neighbours impinge on k
            fpy(list, k) = (flj.*dely); %force * distance
            fpx(k, list) = -fpx(list,k)'; %force that k impinges on neighbours
            fpy(k, list) = -fpy(list,k)'; %transpose column to row
        end
    end
    fx = sum(fpx); %why?
    fy = sum(fpy);
end