close all
clear all
p_gas = [24.3891;46.8891;76.1927;113.172];%gas pressure in mbar
p_air = [18.285;39.9087;69.0084;104.06];%air pressure in mbar
P_therm = [1,1.5,2, 2.5];%thermal power
figure
scatter(p_gas,P_therm)
poly = polyfit(P_therm.',p_gas,2)
pin = polyval(poly,3.1)
y1 = P_therm(1);
m = (P_therm(3) - P_therm(1))/(p_gas(3) - p_gas(1));
x1 = p_gas(1);
Phi_range = 1;
phi = 0.952;%0.952 recommended
Hu_fuel = 9.01*3600;%kW s/m3
p_ch4 = 20:2:160;
P_therm_exp = m.*(p_ch4-x1)+y1; 
vdot_fuel = (P_therm_exp./Hu_fuel);%m3/s
vdot_fuel_lpm = (P_therm_exp./Hu_fuel).*(60000);%lpm
rho_fuel = 0.78;%kg/m3
AF_stoic = 2*(32+3.76*28)/16;
figure
for phi = 0.6:0.05:1
    mdot_air = vdot_fuel.*rho_fuel.*(AF_stoic)/phi;
    vdot_air_lnpm = (mdot_air./1.225).*60000;%lnpm

% for i =1:Phi_range
    plot(p_ch4,vdot_air_lnpm,'LineWidth',2)
    hold on
% end
end
title(sprintf('Air Flow v/s CH4 pressure'));
xlabel('CH4 Pressure(mbar)');
ylabel('Air Flow(lnpm)');
leg=legend(string(0.6:0.05:1), 'Location','northwest');
title(leg,'\phi');
grid on


figure
% for i =1:Phi_range
    plot(p_ch4,P_therm_exp,'LineWidth',2)
    hold on
% end
title(sprintf('Thermal power v/s CH4 pressure'));
xlabel('CH4 Pressure(mbar)');
ylabel('Power(kW)');
grid on

figure
% for i =1:Phi_range
    plot(p_ch4,vdot_fuel_lpm,'LineWidth',2)
    hold on
% end
title(sprintf('CH4 Flow v/s CH4 pressure'));
xlabel('CH4 Pressure(mbar)');
ylabel('CH4 flow(lnpm)');
grid on
