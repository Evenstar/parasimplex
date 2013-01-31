function [soln] = findsoln(n,m,xstar,BB)
   soln = zeros(n,1);
    for i = 1 : m
        if BB(i) <= n
            soln(BB(i)) = xstar(i);
        end
    end
end
