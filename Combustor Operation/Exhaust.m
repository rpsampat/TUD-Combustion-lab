classdef Exhaust<handle
    properties
        T_exh_max = 400+273%K
        Dia = 0.100%m
        rho_air_stp
        mdot_air_cooling = 0.0
        vdot_air_cooling
        T_exh%K
    end
    methods
        function obj = Exhaust(settings,combustor,main_burner)
            obj.rho_air_stp = settings.rho_air_stp;
            obj.temperature(combustor, settings.P_therm,main_burner.mdot_air,main_burner.mdot_fuel)
        end
        function [] =temperature(obj, combustor, P_therm, mdot_air, mdot_fuel)
            cp_air=1;%kJ/kg-K
            cp_fuel = 1.1;
            T_air = 300;%K
            T_fuel = 300;
            cp_mix = 1;
            T_gas = combustor.T_exhaust;
            if obj.mdot_air_cooling == 0.0
%                 obj.T_exh = combustor.T_exhaust;
                if T_gas>obj.T_exh_max
                    P_exh_in = combustor.P_exh;
                    cp_air=1;%kJ/kg-K
                    cp_fuel = 2.3;
                    T_air = 300;%K
                    T_fuel = 300;
                    cp_mix = 1;
                    obj.mdot_air_cooling = (((mdot_air+mdot_fuel)*cp_mix*T_gas) ...
                        -((mdot_air+mdot_fuel)*cp_mix*obj.T_exh_max))/ ...
                        (cp_mix*obj.T_exh_max-cp_air*T_air);
                    %                 obj.mdot_air_cooling = P_exh_in/(cp_air*T_air + cp_air*obj.T_exh_max);
                    obj.vdot_air_cooling = (obj.mdot_air_cooling/obj.rho_air_stp)*60000;%lnpm
                    obj.T_exh = obj.T_exh_max;
                else
                    %                 obj.mdot_air_cooling = 0;
                    obj.vdot_air_cooling = 300;%lnpm
                    obj.mdot_air_cooling = obj.vdot_air_cooling*obj.rho_air_stp/60000;
                    obj.T_exh = (((mdot_air+mdot_fuel)*cp_mix*T_gas) ...
                        +(obj.mdot_air_cooling*cp_air*T_air))/ ...
                ((mdot_air+mdot_fuel+obj.mdot_air_cooling)*cp_mix); 
                end
            else
%                 obj.T_exh = 1000;
                obj.T_exh = (((mdot_air+mdot_fuel)*cp_mix*T_gas) ...
                        +(obj.mdot_air_cooling*cp_air*T_air))/ ...
                ((mdot_air+mdot_fuel+obj.mdot_air_cooling)*cp_mix);                
            end
            
        end
    end
end
