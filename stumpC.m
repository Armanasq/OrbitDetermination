function c = stumpC(z)
    %z - input argument
    %c - value of C(z)

    if z > 0
        c   = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c   = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c   = 1/2;
end
