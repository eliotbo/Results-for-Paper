function M = dctmatrix(d)

for k=0:(d-1)

    for n=0:(d-1)
            if k==0
            alpha = 1/sqrt(2);
            else
                alpha = 1;
             end
        M(k+1,n+1) = sqrt(2/d)*alpha*cos((2*n+1)*k*pi/2/d);
    end
end