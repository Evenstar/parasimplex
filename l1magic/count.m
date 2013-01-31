function out = count(Mp)
[m n] = size(Mp);
grad = Mp(2:end,:)-Mp(1:end-1,:);
s = sign(grad).*(abs(grad)>1e-5);
s2 = s(2:end,:) - s(1:end-1,:);
col = sum(abs(s2),2);
out = sum(sign(col));
end