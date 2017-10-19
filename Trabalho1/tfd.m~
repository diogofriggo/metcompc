function G = tfd(t,f)
    [j,~] = size(f); 
    if(j==1) 
        f=f'; 
    end
    [~,k] = size(t); 
    if(k==1) 
        t=t';
    end
    N=length(t);
    Dt=t(2)-t(1);
    wN=pi/Dt;
    w=linspace(-wN,wN,N)';
    E=w*t;
    %printmat(E)
    W=exp(-1i*E);
    %printmat(imag(W))
    g = W*f;
    printmat(g)
    gn=g/(sqrt(N-1));
    figure(1);plot(t,f);
    figure(2);plot(w,real(gn));
    %figure(3);plot(w,imag(gn))
    %G=[w',gn'];
end