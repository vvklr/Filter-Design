syms s z;                %creating sysmbolic variables
w = linspace(0,2,2001);  %defining Axis Line span from 0 to 2 in 2000steps
jw = linspace(0,2*1i,2001); %defining Axis imaginary 
num = 0.1280;  %first constant of Hanalog LPF
B = 0.1252;
Wo = 0.3199;
N=5;
p = zeros([1 5]);
p(1) = -0.0293 + 0.9926i;
p(2) = -0.2417 + 0.6134i;
p(3) = -0.2987;
p(4) = conj(p(2));
p(5) = conj(p(1));
H_jw = zeros([1 2001]);
H_jw = num + H_jw;
D_s = 1;
for i = 1:5
    H_jw = H_jw./(jw - p(i));
    D_s = D_s*(s - p(i));
end
D_s = sym2poly(D_s)
figure(1);
plot(w,abs(H_jw));


H_jw = num*((Wo^2 - w.^2).^N);
D_s = 1;
N_s = num*((s^2 + Wo^2)^N);
for i = 1:N
    H_jw = H_jw./(p(i)*w.^2 + B.*jw -p(i)*Wo^2);
    D_s = D_s*(-p(i)*s^2 + B*s - p(i)*Wo^2);
end






D_s = sym2poly(D_s)
N_s = sym2poly(N_s)
figure(2);
plot(w, abs(H_jw));



w = linspace(0,pi,1001);
jw = linspace(0,pi*1i,1001);
e_jw = exp(jw);

H_e_jw = num*(((Wo^2 + 1).*((e_jw).^2) + 2.*e_jw.*(Wo^2 - 1) + (Wo^2 + 1)).^N);
D_z = 1;
N_z = num*(((Wo^2 + 1)*(z^2) + 2*z*(Wo^2 - 1) + (Wo^2 + 1))^N);
for i = 1:N
   H_e_jw = H_e_jw./((e_jw.^2).*(B - p(i) - p(i)*Wo^2) - 2.*e_jw.*p(i).*(Wo^2 - 1) - (B + p(i) + p(i)*Wo^2));
   D_z = D_z*((z^2)*(B - p(i) - p(i)*Wo^2) - 2*z*p(i)*(Wo^2 - 1) - (B + p(i) + p(i)*Wo^2));
end


D_z = sym2poly(D_z)
N_z = sym2poly(N_z)
figure(3);
w = w./pi;
plot(w, abs(H_e_jw));