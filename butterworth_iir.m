syms s z;
w = linspace(0,2,2001);
jw = linspace(0,2*1i,2001);
B = 0.1993;
Wo = 0.31023;
Wc = 0.9394;
N = 12;
num= 2.3448;
p = zeros([1 12]);
H_jw = linspace(1,1,2001);
D_s = 1;
N_s = 1;
for i = 1:N
    p(i) = 1i*Wc*exp(1i*(2*i-1)*pi/(2*N));
    H_jw = (H_jw*p(i))./(jw - p(i));
    D_s = D_s*(s - p(i));
    N_s = N_s*p(i);
end
figure(1);
plot(w,abs(H_jw));
N_s;
D_s = sym2poly(D_s);



N_s = num*((B*s)^12);

D_s = 1;
w = linspace(0,5,5001);
jw = linspace(0,5*1i,5001);
H_jw = num*((B.*jw).^12);
for i = 1:N
    H_jw = H_jw./(-w.^2 - p(i)*B.*jw + Wo^2);
    D_s = D_s*(s^2 -p(i)*B*s + Wo^2);
end



D_s = sym2poly(D_s)
N_s = sym2poly(N_s)
figure(2);
plot(w, abs(H_jw));


w = linspace(0,pi,1001);
jw = linspace(0,pi*1i,1001);
e_jw = exp(jw);

H_e_jw = num*((B.*(e_jw.^2 - 1)).^12);
D_z = 1;
N_z = num*B^12*(z^2 - 1)^12;
for i = 1:N
    H_e_jw = H_e_jw./((e_jw.^2).*(1 - p(i)*B + Wo^2) + 2.*e_jw.*(Wo^2 - 1) + (1 + p(i)*B + Wo^2));
    D_z = D_z*((z^2)*(1 - p(i)*B + Wo^2) + 2*z*(Wo^2 - 1) + (1 + p(i)*B + Wo^2));
end
% 
D_z = sym2poly(D_z)
N_z = sym2poly(N_z)
figure(3);
w = w./pi;
plot(w, abs(H_e_jw));