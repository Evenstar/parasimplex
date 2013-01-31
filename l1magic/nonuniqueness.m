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
for i = 0 : 3000
    lambda = 0.01*i + 1e-6;
    tmp = getinstance(obj,lambda);
    Md(i+1,:) = tmp';
    xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-6);
   Mp(i+1,:) = xp';
end

L = sum(abs(Mp),2);
for i = 1 : 4
plot(L, Mp(:,i),'k-');hold all;
end
xlabel('$$\|\beta\|_{l_1}$$','interpreter','latex');
ylabel('Coefficients');
title('primal-dual interior point');
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 pd1.eps
close
count(Mp)
%%
for i = 1 : 2820
    err1(i) = norm(beta-Md(i,:)',1);
    err2(i) = norm(beta-Mp(i,:)',1);
end
L = 1:2820;
L = L * 0.01+ 1e-6;
plot(L,err1,'k-.',L,err2,'k-');
legend('parametric simplex solver','primal-dual interior point');
xlabel('$$\|\beta\|_{l_1}$$','interpreter','latex');
ylabel('$$\|\beta-\hat\beta\|_{l_2}$$','interpreter','latex');
