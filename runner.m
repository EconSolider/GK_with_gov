%% general
clc,clear;
ch=0.42;
pi_ss=1;
z_ss=1;
pistar_ss=1;
d_ss=1;
epsilo=11;
pm_ss=(epsilo-1)/epsilo;
rho_G=0.9;
rho_z=0.9;
rho_i=0.9;
rho_S=0.9;
kapp=0.7;  %GDP成長に対する確率の反応

chi=0.42;  %資本調整コスト
phi_pi=1.5;

delt=0.025;
I_K=delt;

%% default prob
PD_ss=0.005;
S_ss=fsolve(@(S_ss) S_ss/(S_ss+1)-PD_ss,1);
fprintf('PD: %0.4f \n',PD_ss);
fprintf('S: %0.4f \n',S_ss);


%% 金利相関
bet=0.996;
R_ss=1/bet;
QK_ss=1;
G_Y=0.2;
SK_ss=0.005; %金利差(民間) calibrate lambd_K
SB_ss=0.00435; %金利差(国債) calibrate lambd_B
RK_ss=R_ss+SK_ss;
RB_ss=R_ss+SB_ss;
fprintf('RK: %0.4f \n',RK_ss);
fprintf('RB: %0.4f \n',RB_ss);
QB_ss=1;      % calibrate iota in 政府債

%% 実物経済
alph=0.35;
Y_K=(alph*pm_ss)^(-1)*(RK_ss-1+delt);
Ye_K=Y_K;
H_K=(Ye_K)^(1/(1-alph));
w_ss=(1-alph)*pm_ss*Ye_K/H_K;
H_ss=1/3;   %calibrate zeta
K_ss=H_ss/H_K;
Ye_ss=Ye_K*K_ss;
Y_ss=Y_K*K_ss;
G_ss=G_Y*Y_ss;
I_ss=delt*K_ss;
C_ss=Y_ss-I_ss-G_ss;

gamm=0.3;  %消費習慣
lambd_ss=(1-bet*gamm)/(1-gamm)/C_ss;

thet=0.7;
x1_ss=pm_ss*lambd_ss*Y_ss/(1-thet*bet);
x2_ss=lambd_ss*Y_ss/(1-thet*bet);
zet=lambd_ss*w_ss*(1-H_ss); %adjust zeta to match H=1/3;

%% 政府債
sigm=0.05;        %5年国債償還比率
gamm_tau=0.15;     %fiscal rule

Gamm=0.06;        %国債清算余分
iot=(RB_ss/((1-Gamm)*PD_ss+(1-PD_ss))-sigm)/(1-sigm)-1;
fprintf('iotaクーポン額: %0.4f \n',iot);
fprintf('Gamma清算比率: %0.4f \n',Gamm);
B_ss=G_ss/(1+gamm_tau-sigm-(1-sigm)*(iot+1));
fprintf('B/Y: %0.4f \n',B_ss/Y_ss);


%% 純資産相関
L_ss=1.5;             % calibrate eta
n_ss=(K_ss+B_ss)/L_ss;
ps=0.96;           %銀行破産しない確率
et=(1-ps)*n_ss/(K_ss+B_ss); %新たな純資産投入 adjust eta to match L=1.5;

mul_ss=((SK_ss*K_ss+SB_ss*B_ss)/(R_ss*n_ss))/(1+(SK_ss*K_ss+SB_ss*B_ss)/(R_ss*n_ss));
fprintf('流動性制約Largerange multi: %0.4f \n',mul_ss);


A_ss=(1-ps)/(1-mul_ss-ps*R_ss)*R_ss;
Omeg_ss=1-ps+ps*A_ss;
in_ss=R_ss;
D_ss=R_ss*(K_ss+B_ss-n_ss);
T_ss=gamm_tau*B_ss;

lambd_K=1/mul_ss*Omeg_ss*(RK_ss-R_ss);
lambd_B=1/mul_ss*Omeg_ss*(RB_ss-R_ss);
fprintf('A: %0.4f \n',A_ss);
fprintf('Omega: %0.4f \n',Omeg_ss);
fprintf('流動性制約param lambd_K: %0.4f \n',lambd_K);
fprintf('流動性制約param lambd_B: %0.4f \n',lambd_B);


save _param;


%% dyanre program
dynare GK_with_gov














