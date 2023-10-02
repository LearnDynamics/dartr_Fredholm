function y = L_operator(f, tgrid, a, b)

tn = length(tgrid);
y = zeros(tn, 1);
for i = 1:tn
    t = tgrid(i);
    y(i) = integral(@(u) f(u).*exp(-u*t).*u.^(-2), a, b);
end
end