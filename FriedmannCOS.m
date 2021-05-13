clear all;
syms a(t)

t_begin = 0;
t_final = 30;


%Hubble constant at t=0
H0=1;


Omega_r1=0;
Omega_l1=0.735;
Omega_m1=0.265;
Omega_k1=1-Omega_r1-Omega_l1-Omega_m1;

%Initial values
a1_0=1e-100;
a1_dot_0=H0*(Omega_m1/a1_0+Omega_r1/(a1_0^2)+Omega_k1+Omega_l1*(a1_0^2))^(1/2);
%a1_dot_0=1;
inval1=[ a1_0 ; a1_dot_0];

[t1,a1]= ode45(@(t1, a1) fr(t1, a1, Omega_r1, Omega_l1, Omega_m1, Omega_k1, H0), [t_begin t_final], a1_0 );

pl=plot(t1,a1(:,1),'r', 'DisplayName', 'Numerical Solution', 'LineWidth', 2);
hold on;
xlim([0 1])
ylim([0 1])
grid on;

A=load('Ringermacher.txt')
sc=scatter(A(:,2),A(:,3), 100,'k','.', 'DisplayName', 'Ringermacher and Mead Data')

lgd=legend([pl, sc], 'Location', 'northwest')
lgd.FontSize=14;
title('Scale Factor Over Time', 'FontSize', 16)
xlabel('Cosmoligical Time', 'FontSize', 14)
ylabel('a(t)', 'FontSize', 14)
hold off;

function dadt = fr(t, a, Omega_r, Omega_l, Omega_m, Omega_k, H0)
dadt =H0*sqrt(Omega_r/a^2+Omega_m/a+Omega_l*a^2+(Omega_k));
end
