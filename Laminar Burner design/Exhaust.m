classdef Exhaust<handle
    properties
        T_exh_max = 400+273%K
        Dia = 0.100%m
        rho_air_stp
        mdot_air_cooling = 0.0
        vdot_air_cooling
        mdot_gas_tot%kg/s
        T_exh%K
        T_gas%K
        vel_orif%m/s
        dp%pressure drop orifice
    end
    methods
        function obj = Exhaust(settings,combustor,main_burner)
            obj.rho_air_stp = settings.rho_air_stp;
            obj.temperature(combustor, settings.P_therm,main_burner.mdot_air,main_burner.mdot_fuel);
            obj.velocity_orifice(settings.exhaust_orifice_dia);
        end
        function [] =temperature(obj, combustor, P_therm, mdot_air, mdot_fuel)
            cp_air=1;%kJ/kg-K
            cp_fuel = 1.1;
            T_air = 300;%K
            T_fuel = 300;
            cp_mix = 1;
            obj.T_gas = combustor.T_exhaust;
            obj.mdot_gas_tot = mdot_air+mdot_fuel
            if obj.mdot_air_cooling == 0.0
%                 obj.T_exh = combustor.T_exhaust;
                if obj.T_gas>obj.T_exh_max
                    P_exh_in = combustor.P_exh;
                    cp_air=1;%kJ/kg-K
                    cp_fuel = 2.3;
                    T_air = 300;%K
                    T_fuel = 300;
                    cp_mix = 1;
                    obj.mdot_air_cooling = (((mdot_air+mdot_fuel)*cp_mix*obj.T_gas) ...
                        -((mdot_air+mdot_fuel)*cp_mix*obj.T_exh_max))/ ...
                        (cp_mix*obj.T_exh_max-cp_air*T_air);
                    %                 obj.mdot_air_cooling = P_exh_in/(cp_air*T_air + cp_air*obj.T_exh_max);
                    obj.vdot_air_cooling = (obj.mdot_air_cooling/obj.rho_air_stp)*60000;%lnpm
                    obj.T_exh = obj.T_exh_max;
                else
                    %                 obj.mdot_air_cooling = 0;
                    obj.vdot_air_cooling = 300;%lnpm
                    obj.mdot_air_cooling = obj.vdot_air_cooling*obj.rho_air_stp/60000;
                    obj.T_exh = (((mdot_air+mdot_fuel)*cp_mix*obj.T_gas) ...
                        +(obj.mdot_air_cooling*cp_air*T_air))/ ...
                ((mdot_air+mdot_fuel+obj.mdot_air_cooling)*cp_mix); 
                end
            else
%                 obj.T_exh = 1000;
                obj.T_exh = (((mdot_air+mdot_fuel)*cp_mix*obj.T_gas) ...
                        +(obj.mdot_air_cooling*cp_air*T_air))/ ...
                ((mdot_air+mdot_fuel+obj.mdot_air_cooling)*cp_mix);                
            end
            
        end
        
        function [] = velocity_orifice(obj, orifice_dia)
            P_gas = 101325%Pa
            Ru=8.314;
            MW_gas=28.8/1000;
            rho_gas=P_gas*MW_gas/(Ru*obj.T_gas)
            A_orif = 3.14*(orifice_dia^2)/4;
            obj.vel_orif = obj.mdot_gas_tot/(rho_gas*A_orif);   
            beta = orifice_dia/0.115;
            Cd = 0.6;%discharge coefficient of orifice
            obj.dp = ((obj.mdot_gas_tot/(rho_gas*A_orif*Cd))^2)*rho_gas*(1-beta^4)/2;
        end
            
    end
end
