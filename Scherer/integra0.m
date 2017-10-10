dt2 = dt*dt; mdt2 = .5*dt2; h=dt/2;
itmax = ceil(tmax/dt);
forcas0;
for it=1:itmax
    fxi=fx;fyi=fy;
    x = x + vx * dt + mdt2 * fxi;
    y = y + vy * dt + mdt2 * fyi;
    forcas0;
    vx = vx + h * (fx + fxi);
    vy = vy + h * (fy + fyi);
    reflete0;
    plot(x, y, 'o')
    pause(0.01)%
end