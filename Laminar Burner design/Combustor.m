classdef Combustor<handle
    properties
        Dia = 0.155% chamber diameter(m)
        th_q = 0.0035%chamber thickness,(m)
        v_gas%
        T_exhaust%
        T_ad%
        Re%nozzle reynolds number
        P_rad%kw
        P_exh%kW
        mdot_gas%kg/s
        cp_gas=1.3;%kJ/kg-K
        sigma=5.67e-8;%W/m2-K4
        Nu_gas
        Nu_air
        T_quartz_in 
        T_quartz_out
        Q_cond
        T_air_main
    end
    methods
        function obj = Combustor(P_therm, mdot_gas, settings)
            obj.mdot_gas = mdot_gas;
            obj.T_air_main = settings.T_heater;
            obj.flowrates(P_therm, settings);
            
        end
        function [] = flowrates(obj, P_therm, settings)  
            % calculation of temperature based on radiation loss
            A_comb=3.14*obj.Dia*0.49;%m2
            A_x=3.14*(obj.Dia^2)/4;
            epsi=0.3;%0.165;%gas emissivity
            e_Rad=obj.sigma*A_comb*epsi/1000;
            e_exh = (obj.mdot_gas*obj.cp_gas);
            e_rem = (obj.mdot_gas*obj.cp_gas*obj.T_air_main)+P_therm+obj.sigma*A_comb*epsi*300^4/1000;
            p = [e_Rad,0,0,e_exh,-e_rem];
            r = roots(p);
            
            P=1e5;
            Ru=8.314;
            MW_gas=28.8/1000;
            
            %             rho=P*MW_gas/(Ru*t)
            T_samp = imag(r);
            index = find(T_samp==0)
            obj.T_exhaust = max(real(r(index)));%(P_therm/(mdot_gas*cp_gas))+300;
            rho_gas=P*MW_gas/(Ru*obj.T_exhaust);
            obj.T_ad = (P_therm/(obj.mdot_gas*obj.cp_gas))+300;
            obj.P_rad = obj.sigma*A_comb*epsi*(obj.T_exhaust^4)/1000;%kW
            obj.P_exh = (obj.mdot_gas*obj.cp_gas*obj.T_exhaust);%kW
            obj.v_gas = obj.mdot_gas/(rho_gas*3.14*(obj.Dia^2)/4);
            nu = 15.25e-6;
            obj.Re = obj.v_gas*obj.Dia/nu
            Pr = 0.59%from tables
            k = 0.10013;%W/mK, from tables
            obj.Nu_gas =1.67*((obj.Re*Pr/(0.49/0.215))^0.333);
            h_gas = obj.Nu_gas*k/obj.Dia
            Pr_air = 0.69;
            k_air = 0.03128;
            d_coolingchannel =0.015;%m
            Q_cool = 0;%lnpm
            v_cool =(Q_cool/60000)/(3.14*((obj.Dia+2*d_coolingchannel)^2-obj.Dia^2)/4)
            nu_air =1.5e-5;
            if Q_cool ==0 
                % natural convection
                Re_air =40;
                Lc_air =obj.Dia;
                h_air = 1.32*((900-300)/0.49)^0.25
            else
                % forced convection
                Re_air =(v_cool*d_coolingchannel/nu_air);
                Lc_air =d_coolingchannel;%pg 129, Heat and mass transfer databook, C.P.Kothandaraman
                obj.Nu_air = 0.989*(Re_air^0.333)*(Pr_air^0.333);
                h_air = obj.Nu_air*k_air/Lc_air
            end
            
            k_quartz =2.68;%@900C
            R1 = (1/(2*3.14*0.49))*(1/(h_gas*(obj.Dia/2-obj.th_q)));
            R2 = (1/(2*3.14*0.49))*(1*log((obj.Dia/2)/(obj.Dia/2-obj.th_q))/(k_quartz));
            R3 = (1/(2*3.14*0.49))*(1/(h_air*(obj.Dia/2)));
            R = R1+R2+R3;
            Q = ((obj.T_exhaust-(300+10))/R)%W% 
            obj.T_quartz_in = obj.T_exhaust - Q*R1;
            obj.T_quartz_out = obj.T_quartz_in - Q*R2;
            obj.Q_cond = Q/1000;
        end
    end
end