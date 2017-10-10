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