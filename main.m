clear all
close all
clc
%% Calcul numérique
S = 0.0154;
Sn = 5 * 10^-5;
az13 = 0.4753;
az32 = 0.4833;
az20 = 0.9142;
g = 9.81;
Q10 = 3 * 10^-5;
Q20 = 0.5 * 10^-5;
a13 = az13 * Sn * sqrt(2*g);
a32 = az32 * Sn * sqrt(2*g);
a20 = az20 * Sn * sqrt(2*g);
    % Point d'équilibre
H20 = ((Q20 + Q10)/(a20))^2;
H30 = ((Q10)/(a32))^2 + H20;
H10 = ((Q10)/(a13))^2 + H30; 
    % Calcul de Rij
R13 = 2 * ((sqrt(H10-H30))/(a13));
R32 = 2 * ((sqrt(H30-H20))/(a32));
R20 = 2 * ((sqrt(H20))/(a20));

%%
% Partie 4.1
A = [-(1)/(S*R13)   (1)/(S*R13) 0;
      (1)/(S*R13)   -((1)/(S))*((1)/(R13)+(1)/(R32))    (1)/(S*R32);
      0 (1)/(S*R32) -((1)/(S))*((1)/(R32)+(1)/(R20)) ];
B = [1/S    0;
     0   0;
     0   1/S];
%C = [1  0   0;    0   1   0;0   0   1 ];
C = [1  0   0];
D =  0 ;
sys = ss(A,B,C,D);

Vp = eig(A); % sys stable asymptotiquement

Ctrb1e = ctrb(A,B(:,1));
Ctrb2e = ctrb(A,B);
rang1e = rank(Ctrb1e); % Commandable avec 1E
rang2e = rank(Ctrb2e); % Commandable avec 2E

rangInd = rank(Ctrb2e(1:3,1:3)); % la plus grande matrice 3x3 est lineairement independante
                                 % indice E1 = 2 
                                 % indice E2 = 1
% Q4 U1 agit plus sur le sys que U2

%% 
% Partie 4.2
FT = tf(sys);
FTH3 = tf(ss(A,B,[0 1 0],D));
FTH2 = tf(ss(A,B,[0 0 1],D));
% sys stable BIBO
Gain = dcgain(FT);
Pi1 = -0.003718;
Pi2 = -0.2722;
Pi3 = -0.5588;
% subplot(2,1,1)
% step(5e-5*FT(:,1))
% title('H1(t)/U1(t)')
% subplot(2,1,2)
% step(2e-5*FT(:,2))
% title('H1(t)/U2(t)')
% figure(2)
% subplot(2,1,1)
% step(5e-5*FTH3(:,1))
% title('H3(t)/U1(t)')
% subplot(2,1,2)
% step(2e-5*FTH3(:,2))
% title('H3(t)/U2(t)')
% figure(3)
% subplot(2,1,1)
% step(5e-5*FTH2(:,1))
% title('H2(t)/U1(t)')
% subplot(2,1,2)
% step(2e-5*FTH2(:,2))
% title('H2(t)/U2(t)')

% 
% S11 = stepinfo(5e-5*FT(:,1))
% 
% 
% S12 = stepinfo(2e-5*FT(:,2))
% 
% 
% S31 = stepinfo(5e-5*FTH3(:,1))
% 
% S32 = stepinfo(2e-5*FTH3(:,1))
% 
% 
% S21 = stepinfo(5e-5*FTH2(:,2))
% 
% 
% S22 = stepinfo(2e-5*FTH2(:,1))

%% 
% Partie 5.2
Vpdes = [-3.5/90 -15/90 -17/90];
q = [1;1];
Btild = B*q;
ctrbBass = ctrb(A,Btild);
rangBass = rank(ctrbBass);
k = acker(A,Btild,Vpdes);
K = q*k ;
sysBass = ss(A-B*K,B,C,D);
gainSBass = dcgain(sysBass);
N = gainSBass'*inv(gainSBass*gainSBass');
S = stepinfo(0.05*sysBass);
tpsReponse = S.SettlingTime
%%
% Partie 5.3
MatCo = ctrb(A,B);
Cob = ComputeCob(A,B);
P = FindMatriceP(A,B);
[Ac,Bc] = ControllableFormMIMO(A,B);
R = FindMatricePre(A,B);
Bcc = Bc * R;
Kco = place(Ac,Bcc,Vpdes)
%Kco = [0.00648 0.1692 -0.0397;0.0001 0 0.1385];
Vpco = eig(Ac-Bcc*Kco);
sysCo = ss(Ac-Bcc*Kco,Bcc,C,D);
gainCo = dcgain(sysCo);
Nco = gainCo'*inv(gainCo*gainCo');
Sco = stepinfo(0.05*sysCo);
tpsReponseCo = Sco.SettlingTime
tpsMonteCo = Sco.RiseTime