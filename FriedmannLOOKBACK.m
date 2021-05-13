clear all;
syms a(t)

%Need forward and backward times
t_begin(1) = 0;
t_final(1) = -30;
t_begin(2) = 0;
t_final(2) = 30;

%Hubble constant today
H0=71000/(3*10^(22))*(3600*24*365*10^(9));

%Cosmological Parameters
Omega_r1=8.4*10^(-4);
Omega_l1=0.735;
Omega_m1=0.265;
Omega_k1=1-Omega_r1-Omega_l1-Omega_m1;

%Initial values
a1_0=1;
a1_dot_0=H0*(1-Omega_k1+Omega_m1+Omega_r1+Omega_l1)^(1/2);
inval1=[ a1_0 ; a1_dot_0];
%Solving forward and back (today forward, today backward)
for i=1:2;
    [t1,a1]= ode45(@(t1, a1) fr(t1, a1, Omega_r1, Omega_l1, Omega_m1, Omega_k1, H0), [t_begin(i) t_final(i)], a1_0 );

    pl(i)=plot(t1,a1(:,1),'r', 'DisplayName', 'Numerical Solution', 'LineWidth', 2);
    hold on;
    xlim([-14 0.2])
    ylim([0 1.2])
    grid on;
end;
A=load('Ringermacher.txt')
%Ringermacher data is in terms of Cosmological Time
%so it must be rescaled for actual time (Universe age = 13.8 Gyrs)
sc=scatter(-13.8+13.8*A(:,2),A(:,3), 100,'k','.', 'DisplayName', 'Ringermacher and Mead Data')

lgd=legend([pl(1), sc], 'Location', 'northwest')
lgd.FontSize=14;
title('Scale Factor Over Time', 'FontSize', 16)
xlabel('t, Lookback Time (Gyr)', 'FontSize', 14)
ylabel('a(t)', 'FontSize', 14)
hold off;

function dadt = fr(t, a, Omega_r, Omega_l, Omega_m, Omega_k, H0)
dadt =H0*sqrt(Omega_r/a^2+Omega_m/a+Omega_l*a^2+(Omega_k));
end






