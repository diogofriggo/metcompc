function dinamol0 (X, Y, l, N, dt, tmax, T)
    %help dinamol0; clear;
    %inicia0;
    while(length(tmax) > 0)
        %integra0;
        Ec = (sum(vx.^2)+sum(vy))
        %potencial0;
        Energia = V+Ec;
        EcEpE = [Ec, V, Energia]
        T = Ec/N
        tmax = input('novo tmax? se "enter", encerra: \n');
        if(length(tmax) > 0);
            renorm = sqrt(Tnova/T);
            vx = vx*renorm;
            vx = vy*renorm;
            T = Tnova
        end
    end
end