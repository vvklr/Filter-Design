wc1 = 0.185*pi;
wc2 = 0.2430*pi;
N1 = 103;
n1 = linspace(-N1, N1, 2*N1 + 1);
h1_n = (sin(wc2.*n1) - sin(wc1.*n1))./(pi.*n1); % defining vector
h1_n(N1 + 1) = (wc2 - wc1)/pi;  %defining Zero at N+1
w = linspace(0, pi, 1001);
jw = linspace(0, pi*1i, 1001);
e_jw = exp(jw);
H1_e_jw = zeros([1 1001]);
for i = 1:(2*N1 + 1)
    H1_e_jw = H1_e_jw + h1_n(i).*e_jw.^(N1 + 1 -i);%z transform
end
figure(1);
plot(w/pi, abs(H1_e_jw));