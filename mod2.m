function y = mod2(a,b)
    y = mod(a,b);
    if y == 0
        y = b;
    end
end