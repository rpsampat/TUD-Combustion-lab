close all
clear all
% input settings
type = "Main burner";%Non-reacting, Pilot burner, Main burner
equivalence_ratio = 0.625;
power_fuel = 60;%kW
P_min = 15;
P_max = 200;
p_step = 5;
P_Range = length(P_min:p_step:P_max);
eq_min = 0.6;
eq_max = 0.85;
eq_step = 0.05;          
Phi_range = length(eq_min:eq_step:eq_max);
P_therm = zeros(Phi_range,P_Range);
P_rad = zeros(Phi_range,P_Range);
P_cond = zeros(Phi_range,P_Range);
T_gas = zeros(Phi_range,P_Range);
T_gas_ad = zeros(Phi_range,P_Range);
T_quartz = zeros(Phi_range,P_Range);
vdot_air = zeros(Phi_range,P_Range);
vdot_air_cool = zeros(Phi_range,P_Range);
v_nozzle = zeros(Phi_range,P_Range);
v_nozzle_fuel = zeros(Phi_range,P_Range);
volflow_fuel = zeros(Phi_range,P_Range);
v_gas = zeros(Phi_range,P_Range);
T_heater = zeros(Phi_range,P_Range);
T_exhaust = zeros(Phi_range,P_Range);
num_equi = 0;
O2_perc = 20.95;

%% varying equivalence ratio(legend) and power(x axis). Dilution =0.
for equivalence_ratio = eq_min:eq_step:eq_max
    num_equi= num_equi+1;
    num_power = 0;
    for power_fuel =P_min:p_step:P_max
        num_power = num_power + 1;
        settings = Settings(power_fuel, equivalence_ratio, O2_perc);
        main_burner = MainBurner(settings)
        mdot_gas = main_burner.mdot_air + main_burner.mdot_fuel;
        comb = Combustor(main_burner.P_therm, mdot_gas, settings);
        exhaust = Exhaust(settings,comb, main_burner);
        operation = Operation(settings, main_burner, comb, exhaust)  ;
        airline = AirLine(settings,operation);
        exhaust.mdot_air_cooling = airline.vdot_FCV001*settings.rho_air_stp/60000;%kg/s
        exhaust.temperature(comb, settings.P_therm, main_burner.mdot_air, main_burner.mdot_fuel)
        heater = Heater(main_burner.mdot_air);
        
        P_therm(num_equi,num_power) = main_burner.P_therm;
        P_rad(num_equi,num_power) = comb.P_rad;
        P_cond(num_equi,num_power) = comb.Q_cond;
        T_gas(num_equi,num_power) = comb.T_exhaust;
        T_gas_ad(num_equi,num_power) = comb.T_ad;
        T_quartz(num_equi,num_power) = comb.T_quartz_in;
        vdot_air(num_equi,num_power) = airline.vdot_FCV002;
        vdot_air_cool(num_equi,num_power) = airline.vdot_FCV001;
        v_nozzle(num_equi,num_power) = main_burner.v_nozzle;
        v_nozzle_fuel(num_equi,num_power) = main_burner.v_nozzle_fuel;
        volflow_fuel(num_equi,num_power) = main_burner.vdot_fuel;
        v_gas(num_equi,num_power) = comb.v_gas;
        T_heater(num_equi,num_power) = heater.Temp;        
        T_exhaust(num_equi,num_power) = exhaust.T_exh;
    end
end

dia_title="6";
figure
for i =1:Phi_range
    plot(P_min:p_step:P_max,P_therm(i,:),'LineWidth',2)
    hold on    
end
title(sprintf('Thermal Power v/s Power Demand'));
xlabel('Power(kW)');
ylabel('Thermal Power(kW)');  
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on


figure
for i =1:Phi_range
    plot(P_therm(i,:),P_rad(i,:),'LineWidth',2)
    hold on    
end
title(sprintf('Radiated Power v/s Power'));
xlabel('Power(kW)');
ylabel('Radiated Power(kW)');  
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on 
saveas(gcf,'radiatedpower_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),P_cond(i,:),'LineWidth',2)
    hold on    
end
title(sprintf('Conducted Power v/s Power'));
xlabel('Power(kW)');
ylabel('Conducted Power(kW)');  
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on 
saveas(gcf,'conductedpower_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),100.*(P_cond(i,:)+P_rad(i,:))./P_therm(i,:),'LineWidth',2)
    hold on    
end
title(sprintf('Heat loss v/s Power'));
xlabel('Power(kW)');
ylabel('Heat loss(%)');  
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northeast');
title(leg,'\phi');
grid on 
saveas(gcf,'heatloss_'+dia_title+'mmnozzle.png')

figure
% plot(P_therm,T_gas,'LineWidth',2)
% hold on
% plot(P_therm,T_gas_ad,'LineWidth',2)
% hold on
for i =1:Phi_range
    plot(P_therm(i,:),T_gas(i,:),'LineWidth',2)
    hold on    
end
% for i =1:Phi_range
%     plot(P_therm(i,:),T_gas_ad(i,:),'LineWidth',2,'LineStyle','--')
%     hold on    
% end
title(sprintf('Gas Temperature v/s Power'));
xlabel('Power(kW)');
ylabel('Temperature(K)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'gastemp_'+dia_title+'mmnozzle.png')

figure
% plot(P_therm,T_gas,'LineWidth',2)
% hold on
% plot(P_therm,T_gas_ad,'LineWidth',2)
% hold on
for i =1:Phi_range
    plot(P_therm(i,:),T_exhaust(i,:),'LineWidth',2)
    hold on    
end
% for i =1:Phi_range
%     plot(P_therm(i,:),T_gas_ad(i,:),'LineWidth',2,'LineStyle','--')
%     hold on    
% end
title(sprintf('Exhaust Temperature v/s Power'));
xlabel('Power(kW)');
ylabel('Temperature(K)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on


figure
% plot(P_therm,T_gas,'LineWidth',2)
% hold on
% plot(P_therm,T_gas_ad,'LineWidth',2)
% hold on
for i =1:Phi_range
    plot(P_therm(i,:),T_quartz(i,:),'LineWidth',2)
    hold on    
end
plot(P_therm(i,:),(950+273)*ones(P_Range,1),'LineWidth',2,'Color','black','LineStyle','--')
hold on
title(sprintf('Quartz Wall Temperature v/s Power'));
xlabel('Power(kW)');
ylabel('Temperature(K)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'quartztemp_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),vdot_air(i,:),'LineWidth',2)
    hold on
end
title(sprintf('Main Air flow rate v/s Power'));
xlabel('Power(kW)');
ylabel('Flow rate(lnpm)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'mainairflow_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),vdot_air_cool(i,:),'LineWidth',2)
    hold on
end
% plot(P_therm,(300)*ones(length(P_therm),1),'LineWidth',2)
% hold on
title(sprintf('Cooling Air flow rate v/s Power'));
xlabel('Power(kW)');
ylabel('Flow rate(lnpm)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on

figure
for i =1:Phi_range
    plot(P_therm(i,:),volflow_fuel(i,:),'LineWidth',2)
    hold on
end
% plot(P_therm,(300)*ones(length(P_therm),1),'LineWidth',2)
% hold on
title(sprintf('Fuel flow rate v/s Power'));
xlabel('Power(kW)');
ylabel('Flow rate(lnpm)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'Fuelflow_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),volflow_fuel(i,:)./settings.nozzles,'LineWidth',2)
    hold on
end
% plot(P_therm,(300)*ones(length(P_therm),1),'LineWidth',2)
% hold on
title(sprintf('Fuel flow rate per nozzle v/s Power'));
xlabel('Power(kW)');
ylabel('Flow rate(lnpm)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
%saveas(gcf,'Fuelflow_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),v_nozzle(i,:),'LineWidth',2)
    hold on
end
title(sprintf('Nozzle Velocity v/s Power'));
xlabel('Power(kW)');
ylabel('Velocity(m/s)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'nozzlevelocity_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),v_nozzle(i,:).*settings.nozzle_dia./settings.nu,'LineWidth',2)
    hold on
end
title(sprintf('Re v/s Power'));
xlabel('Power(kW)');
ylabel('Re');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'Re_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),v_nozzle_fuel(i,:),'LineWidth',2)
    hold on
end
title(sprintf('Fuel Nozzle Velocity v/s Power'));
xlabel('Power(kW)');
ylabel('Velocity(m/s)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'fuelnozzlevelocity_'+dia_title+'mmnozzle.png')


figure
for i =1:Phi_range
    plot(P_therm(i,:),v_gas(i,:),'LineWidth',2)
    hold on
end
title(sprintf('Combustor Velocity v/s Power'));
xlabel('Power(kW)');
ylabel('Velocity(m/s)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'combustorvelocity_'+dia_title+'mmnozzle.png')

figure
for i =1:Phi_range
    plot(P_therm(i,:),T_heater(i,:)-273,'LineWidth',2)
    hold on
end
title(sprintf('Heater temperature v/s Power'));
xlabel('Power(kW)');
ylabel('Temperature(C)');
leg=legend(string(eq_min:eq_step:eq_max), 'Location','northwest');
title(leg,'\phi');
grid on
saveas(gcf,'heatertemp_'+"48kW_"+dia_title+'mmnozzle.png')



    

