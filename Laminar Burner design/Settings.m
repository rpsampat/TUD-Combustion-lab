classdef Settings<handle
    properties
        nozzles% number of nozzles
        P_therm%thermal power(kW)
        LHV%Lower heating value of fuel(MJ/kg)
        nozzle_dia%m
        nozzle_dia_fuel%m
        phi_main%Main burner equivalence ratio
        MW_fuel = 16/1000%kg/mol
        AF_stoic = 2*(32+3.76*28)/16;%for CH4/air
        MW_air = 28.8/1000%kg/mol
        MW_N2 = 28/1000%kg/mol
        nu = 1.5e-5
        rho_air_stp
        P_CAV001 = 40e5;%Air tank pressure, Pa
        T_CAV001 = 273.15+ 15;%Air tank temperature. K   
        T_heater = 273.15+ 15;
        O2_perc
        exhaust_orifice_dia
    end
    methods
        function obj = Settings(power_fuel, phi, O2_perc)
            obj.nozzles = 1;% number of nozzles
            obj.P_therm = power_fuel;
            obj.LHV = 50;%MJ/kg
            obj.nozzle_dia = 0.021;%m
            obj.nozzle_dia_fuel = 0.001;%m
            obj.phi_main = phi;
            obj.rho_air_stp = 101300*obj.MW_air/(8.314*273.15);%kg/m3   
            obj.O2_perc = O2_perc;%percentage O2 by volume
            obj.exhaust_orifice_dia = 0.040;%m
        end
    end
end
