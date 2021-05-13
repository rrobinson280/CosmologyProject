clear all;
syms a(t)

t_begin(1) = 0;
t_final(1) = 4;


%Hubble's constant at t=0
H0=1;

%Cosmological Parameters
Omega_r1=0;
Omega_l1=0;
Omega_m1=1;
Omega_k1=1-Omega_r1-Omega_l1-Omega_m1;

Omega_r2=0;
Omega_l2=0;
Omega_m2=5;
Omega_k2=1-Omega_r2-Omega_l2-Omega_m2;

Omega_r3=8e-4;
Omega_l3=0.7;
Omega_m3=0.3;
Omega_k3=1-Omega_r3-Omega_l3-Omega_m3;

Omega_r4=0;
Omega_l4=0;
Omega_m4=-2;
Omega_k4=1-Omega_r4-Omega_l4-Omega_m4;

%Initial values
a1_0=1e-100;
a1_dot_0=sqrt(H0*(Omega_m1/a1_0+Omega_r1/(a1_0^2)+Omega_k1+Omega_l1*(a1_0^2))^(1/2));
inval1=[ a1_0 ; a1_dot_0];

%Solving the different universes
[t1,a1]= ode45(@(t1, a1) fr(t1, a1, Omega_r1, Omega_l1, Omega_m1, Omega_k1, H0), [t_begin t_final], [a1_0; a1_dot_0]);
[t2,a2]= ode45(@(t2, a2) fr(t2, a2, Omega_r2, Omega_l2, Omega_m2, Omega_k2, H0), [t_begin 0.9445], [a1_0; a1_dot_0]);
[t3,a3]= ode45(@(t3, a3) fr(t3, a3, Omega_r3, Omega_l3, Omega_m3, Omega_k3, H0), [t_begin t_final], [a1_0; a1_dot_0] );
[t4,a4]= ode45(@(t4, a4) fr(t4, a4, Omega_r4, Omega_l4, Omega_m4, Omega_k4, H0), [t_begin t_final], [a1_0; a1_dot_0] );

%Plotting
plm=plot(t1,a1(:,1),'g', 'DisplayName', 'k=0, Einstein-de Sitter', 'LineWidth', 2);
hold on;
pl=plot(t2,a2(:,1),'r', 'DisplayName', 'k=1', 'LineWidth', 2);
pl2=plot(-(t2-2*0.9445),a2(:,1),'r', 'LineWidth', 2);
pl3=plot(t3,a3(:,1), 'b', 'DisplayName', 'k=0, \Omega_{0m}=0.4, \Omega_{0r}=8e-4, \Omega_{0\Lambda}=0.7', 'LineWidth', 2);
pl4=plot(t4,a4(:,1), 'c', 'DisplayName','k=-1', 'LineWidth', 2); 
xlim([0 4]);
ylim([0 4]);
grid on;

lgd=legend([pl, plm, pl3, pl4], 'Location', 'southeast')
lgd.FontSize=11;
title('Scale Factor Over Time', 'FontSize', 16)
xlabel('Cosmoligical Time', 'FontSize', 14)
ylabel('a(t)', 'FontSize', 14)
hold off;

function dadt = fr(t, a, Omega_r, Omega_l, Omega_m, Omega_k, H0)
dadt = zeros(2,1);
dadt(1) = H0*(Omega_r*a(1)^(-2)+Omega_m*a(1)^(-1)+Omega_l*a(1)^2+(Omega_k))^(1/2);
dadt(2) = H0^2/2*(-Omega_m*a(1)^(-2)-Omega_r*a(1)^(-3)+2*Omega_l*a(1));
end





