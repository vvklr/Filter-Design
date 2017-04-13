%BANDpass
omega1=11.5*2*pi/140;
omega2=16.5*2*pi/140;
T_wid=2*pi/140;
B=omega2-omega1;
omega0=sqrt(omega1*omega2);
a=((omega1-T_wid)^2-(omega0)^2)/(B*(omega1-T_wid));
b=((omega2+T_wid)^2-(omega0)^2)/(B*(omega2+T_wid));
omegas=min(abs(a),b);
del2=0.1;
del1=0.1;
D2=1/(del2^2)-1;
D1=(1/((1-del1)^2))-1;
N=ceil(acosh(sqrt(D2)/sqrt(D1))/acosh(omegas/1));
% omegac=1*cosh((1/N)*(acosh(sqrt((1+2*D1)/D1))));
den=cosh(N*acosh(omegas));
epsilonl=sqrt(D2/den);
epsilonu=sqrt(D1);

epsilon = (epsilonl+epsilonu)/2;

 H_bp = tf([1]/[1]);
 H_lp = tf([1]/[1]);
     for k = 1:N
          % pole= w_c*exp((1i*(2*k+N-1)*pi)/(2*N));
         pole=-sinh((1/N)*asinh(1/epsilon))*sin((pi/2)*(2*k-1)/N)+1j*cosh((1/N)*asinh(1/epsilon))*cos((pi/2)*(2*k-1)/N);
         
         H_bp = H_bp* tf([B*pole 0],[1 -pole*B +(omega0^2)]);
         H_lp = H_lp*tf([1],[1 -pole]);
     end
 H_lp = tf(real(cell2mat(H_lp.num)), real(cell2mat(H_lp.den)));
 H_bp = tf(real(cell2mat(H_bp.num)), real(cell2mat(H_bp.den)));
 [num, den] = bilinear(cell2mat(H_bp.num), cell2mat(H_bp.den),0.5);

 freqz(num,den);

% %%%%%%%%%%%BANDSTOP IIR FILTER DESIGN CHEBYSHEV%%%%%%%%%%%
% 
% omega1s=9.5*2*pi/140;
% omega2s=12.5*2*pi/140;
% Bs=omega2s-omega1s;
% omega0s=sqrt(omega1s*omega2s);
% as=(Bs*omega1s+T_wid)/((-(omega1s+T_wid)^2+(omega0s)^2));
% bs=(Bs*omega2s-T_wid)/((-(omega2s-T_wid)^2+(omega0s)^2));
% 
% omegasb=min(abs(as),abs(bs));
% 
% 
% Nbs=ceil(0.5*log(D2/D1)/log(omegasb/1));
% omegacl=1*((D1)^(-0.5*Nbs));
% omegacb=omegasb*(D2^(-0.5*Nbs));