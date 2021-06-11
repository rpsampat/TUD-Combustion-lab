classdef MainBurner<handle
    properties
        N% number of nozzles
        P_therm%thermal power(kW)
        LHV%Lower heating value of fuel(MJ/kg)
        phi%equivalence ratio
        vdot_air%air flow lnpm
        vdot_fuel%fuel flow lnpm
        vdot_fuel_min = 0.01%minimum fuel flow lnpm
        vdot_fuel_max = 350%maximum fuel flow lnpm
        nozzle_dia%m
        nozzle_dia_fuel%m
        v_nozzle%velocity at nozzle m/s
        v_nozzle_fuel
        Re%nozzle reynolds number
        mdot_fuel% mass flow fuel kg/s
        mdot_air% mass flow air kg/s
        T_air_main
        O2_perc %percentage O2 by volume
        mdot_N2
        vdot_N2
    end
    methods
        function obj = MainBurner(settings)
            obj.N = settings.nozzles;
            obj.P_therm = settings.P_therm;
            obj.phi = settings.phi_main;
            obj.LHV = settings.LHV;
            obj.nozzle_dia = settings.nozzle_dia;
            obj.nozzle_dia_fuel = settings.nozzle_dia_fuel;
            obj.T_air_main = settings.T_heater;
            
            obj.flowrates(settings);
            
        end
        function [] = flowrates(obj, settings)
            mdot_fuel = (obj.P_therm./obj.LHV)*(1/1000);%kg/s
            rho_fuel_stp = (101300)*settings.MW_fuel/(8.314*273.15);%kg/m3
            obj.vdot_fuel = (mdot_fuel/rho_fuel_stp)*(60*1000);%lnpm    
            if obj.vdot_fuel< obj.vdot_fuel_min
                obj.vdot_fuel = obj.vdot_fuel_min;%lnpm
                obj.P_therm = (obj.vdot_fuel/60000)*rho_fuel_stp*obj.LHV*1000;%kW
            elseif obj.vdot_fuel> obj.vdot_fuel_max
                obj.vdot_fuel = obj.vdot_fuel_max;%lnpm
                obj.P_therm = (obj.vdot_fuel/60000)*rho_fuel_stp*obj.LHV*1000;%kW
            end
            obj.mdot_fuel = obj.vdot_fuel*rho_fuel_stp/60000;
            obj.vdot_air = (obj.vdot_fuel*settings.AF_stoic/obj.phi)*(settings.MW_fuel/settings.MW_air);%lnpm
            rho_air_stp = 101300*settings.MW_air/(8.314*273.15);%kg/m3
            rho_air_main = 101300*settings.MW_air/(8.314*obj.T_air_main);%kg/m3
            rho_N2_stp = 101300*settings.MW_N2/(8.314*273.15);%kg/m3
            obj.mdot_air = obj.vdot_air * rho_air_stp/60000;
            obj.mdot_N2 = (20.95*obj.mdot_air/settings.O2_perc)-obj.mdot_air;
            obj.vdot_N2 = obj.mdot_N2*60000/rho_N2_stp;%lnpm
            mdot_nozzle = (obj.mdot_air+obj.mdot_N2+obj.mdot_fuel)/obj.N;
            mdot_nozzle_fuel = (obj.vdot_fuel*rho_fuel_stp/60000)/obj.N;
            area_nozzle = 3.14*(obj.nozzle_dia)^2/4;
            area_nozzle_fuel = 3.14*(obj.nozzle_dia_fuel)^2/4;
            obj.v_nozzle = (mdot_nozzle)/(rho_air_main*area_nozzle);%m/s
            obj.v_nozzle_fuel = mdot_nozzle_fuel/(rho_fuel_stp*area_nozzle_fuel);%m/s
            obj.Re = obj.v_nozzle*obj.nozzle_dia/settings.nu;
            
        end
    end
end

        
            
            
           