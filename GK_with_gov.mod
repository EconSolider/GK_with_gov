var
C lambd H w R
Omeg A RK RB mul
n QK QB K B
D L SK SB Ye
pm pistar x1 x2 z
pi PD S T G
I in Y d
;

varexo 
 eps_z eps_i eps_G eps_S
;

parameters
ch gamm zet bet ps et alph delt epsilo thet
kapp Gamm gamm_tau sigm rho_G chi phi_pi rho_i rho_z rho_S 
lambd_K lambd_B iot
%steady_state_values
C_ss lambd_ss H_ss w_ss R_ss
Omeg_ss A_ss RK_ss RB_ss mul_ss
n_ss QK_ss QB_ss K_ss B_ss
D_ss L_ss SK_ss SB_ss Ye_ss
pm_ss pistar_ss x1_ss x2_ss z_ss
pi_ss PD_ss S_ss T_ss G_ss
I_ss in_ss Y_ss d_ss
;

load '_param.mat';
set_param_value('ch',ch);
set_param_value('gamm',gamm);
set_param_value('zet',zet);
set_param_value('bet',bet);
set_param_value('ps',ps);
set_param_value('et',et);
set_param_value('alph',alph);
set_param_value('delt',delt);
set_param_value('epsilo',epsilo);
set_param_value('thet',thet);
set_param_value('kapp',kapp);
set_param_value('Gamm',Gamm);
set_param_value('gamm_tau',gamm_tau);
set_param_value('sigm',sigm);
set_param_value('rho_G',rho_G);
set_param_value('chi',chi);
set_param_value('phi_pi',phi_pi);
set_param_value('rho_i',rho_i);
set_param_value('rho_z',rho_z);
set_param_value('rho_S',rho_S);
set_param_value('iot',iot);
set_param_value('lambd_K',lambd_K);
set_param_value('lambd_B',lambd_B);

set_param_value('C_ss',C_ss);
set_param_value('lambd_ss',lambd_ss);
set_param_value('H_ss',H_ss);
set_param_value('w_ss',w_ss);
set_param_value('R_ss',R_ss);
set_param_value('Omeg_ss',Omeg_ss);
set_param_value('A_ss',A_ss);
set_param_value('RK_ss',RK_ss);
set_param_value('RB_ss',RB_ss);
set_param_value('mul_ss',mul_ss);
set_param_value('n_ss',n_ss);
set_param_value('QK_ss',QK_ss);
set_param_value('QB_ss',QB_ss);
set_param_value('K_ss',K_ss);
set_param_value('B_ss',B_ss);
set_param_value('D_ss',D_ss);
set_param_value('L_ss',L_ss);
set_param_value('SK_ss',SK_ss);
set_param_value('SB_ss',SB_ss);
set_param_value('Ye_ss',Ye_ss);
set_param_value('pm_ss',pm_ss);
set_param_value('pistar_ss',pistar_ss);
set_param_value('x1_ss',x1_ss);
set_param_value('x2_ss',x2_ss);
set_param_value('z_ss',z_ss);
set_param_value('pi_ss',pi_ss);
set_param_value('PD_ss',PD_ss);
set_param_value('S_ss',S_ss);
set_param_value('T_ss',T_ss);
set_param_value('G_ss',G_ss);
set_param_value('I_ss',I_ss);
set_param_value('in_ss',in_ss);
set_param_value('Y_ss',Y_ss);
set_param_value('d_ss',d_ss);


model;
1/(C-gamm*C(-1))-lambd-bet*gamm/(C(+1)-gamm*C)=0;
lambd*w=zet/(1-H);
bet*lambd(+1)-lambd/R=0;
Omeg=lambd(+1)/lambd*(1-ps+ps*A(+1)); %%
Omeg*(RK(+1)-R)=mul*lambd_K; %%
Omeg*(RB-R)=mul*lambd_B; %%
A=Omeg*R/(1-mul); %%
mul=1-Omeg*R*n/(lambd_K*QK*K+lambd_B*QB*B);
D/R=QK*K+QB*B-n; %%
L=(QK*K+QB*B)/n;
SK=RK(+1)-R;
SB=RB-R;
n=ps*n(-1)+et*(QK(-1)*K(-1)+QB(-1)*B(-1));
Ye=z*K(-1)^alph*H^(1-alph);
w*H=(1-alph)*pm*Ye;
RK=alph*pm*Ye/(QK(-1)*K(-1))+QK*(1-delt)/QK(-1);
pistar=epsilo/(epsilo-1)*x1/x2;
x1=lambd*pm*Y+thet*bet*pi(+1)^epsilo*x1(+1);
x2=lambd*Y+thet*bet*pi(+1)^(epsilo-1)*x2(+1);
1=thet*pi^(epsilo-1)+(1-thet)*pistar^(1-epsilo);
PD=S/(1+S);   %%
log(S)=(1-rho_S)*log(steady_state(S))+rho_S*log(S(-1))+eps_S;
RB=((1-Gamm)*PD(+1)+(1-PD(+1)))*(sigm+(1-sigm)*(iot+QB(+1)))/QB; %%
T=gamm_tau*B(-1);
QB*B=(sigm+(1-sigm)*(iot+QB))*B(-1)+G-T;
log(G)=rho_G*log(G(-1))+(1-rho_G)*log(steady_state(G))+eps_G;
K=(1-delt)*K(-1)+I-ch/2*((I/K(-1)-delt)^2)*K(-1);
QK=(1-ch*(I/K(-1)-delt))^(-1);
in=R/pi;
log(in/steady_state(in))=phi_pi*log(pi/steady_state(pi))+eps_i;
Ye=d*Y;
d=thet*(pi^epsilo)*d(-1)+(1-thet)*pistar^(-epsilo);
log(z)=rho_z*log(z(-1))+eps_z;
Y=C+I+G+ch/2*((I/K(-1)-delt)^2)*K(-1);
end;

steady_state_model;
C=C_ss;
lambd=lambd_ss;
H=H_ss;
w=w_ss;
R=R_ss;
Omeg=Omeg_ss;
A=A_ss;
RK=RK_ss;
RB=RB_ss;
mul=mul_ss;
n=n_ss;
QK=QK_ss;
QB=QB_ss;
K=K_ss;
B=B_ss;
D=D_ss;
L=L_ss;
SK=SK_ss;
SB=SB_ss;
Ye=Ye_ss;
pm=pm_ss;
pistar=pistar_ss;
x1=x1_ss;
x2=x2_ss;
z=z_ss;
pi=pi_ss;
PD=PD_ss;
S=S_ss;
T=T_ss;
G=G_ss;
I=I_ss;
in=in_ss;
Y=Y_ss;
d=d_ss;
end;

steady;
check;

shocks;
var eps_z=0.01;
var eps_S=0.01;
var eps_i=0.01;
end;

stoch_simul(order=1,irf=20,periods=0);




























