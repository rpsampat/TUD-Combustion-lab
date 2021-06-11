classdef Heater<handle
    properties
        P_heater = 48;%kW
        m_air
        Cp = 1.002;%kJ/kg-K
        Temp
    end
    methods
        function obj = Heater(m_air)
            obj.m_air = m_air;
            obj.Temperature();
        end
        function [] = Temperature(obj)
            obj.Temp = obj.P_heater/(obj.m_air*obj.Cp);
            if obj.Temp>(700+273)
                obj.Temp = 700+273;
            end
        end
    end
end