%%
clear all
X = [1 1 2 1;
     1 0 1 2;
     0 1 1 0];
beta = [2 0 3 1]';
Y = [9.1 7.2 2.8]';
obj = dantzig(X,Y,1e-6);
figure(1);
L = sum(abs(obj.M),2);
for i = 1 : 4
plot(L, obj.M(:,i),'k-');hold all;
end
xlabel('$$\|\beta\|_{l_1}$$','interpreter','latex');
ylabel('Coefficients');
title('parametric simplex solver');
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 ps1.eps
close;
figure;
x0 = X'*Y;
for i = 0 : 300
    lambda = 0.1*i + 1e-6;
    xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-4);
   Mp(i+1,:) = xp';
end

L = sum(abs(Mp),2);
for i = 1 : 4
plot(L, Mp(:,i),'k-');hold all;
end
xlabel('$$\|\beta\|_{l_1}$$','interpreter','latex');
ylabel('Coefficients');
title('l1 magic');
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 magic1.eps
close
%%
