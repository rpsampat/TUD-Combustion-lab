classdef AirLine<handle
    properties
        P_CAV001%Air tank pressure
        T_CAV001%Air tank temperature
        vdot_LN001%Total air flow
        vdot_PCV001%high flow pressure reducer
        vdot_PCV002%low flow pressure reducer
        PCV002_lowlim = 80;%lnpm
        PCV002_highlim = 3400;%lnpm
        PCV001_lowlim = 9700;%lnpm
        PCV001_highlim = 30000;%lnpm
        P_PT003
        vdot_FCV001%cooling air
        vdot_FCV001_lowlim = 300%lnpm
        vdot_FCV001_highlim = 15000%lnpm
        vdot_FCV002%main burner air
        vdot_FCV002_lowlim = 300%lnpm
        vdot_FCV002_highlim = 15000%lnpm
        vdot_FCV005%pilot burner air
    end
    methods
        function obj = AirLine(settings,operation)
            obj.P_CAV001 = settings.P_CAV001;%Air tank pressure
            obj.T_CAV001 = settings.T_CAV001;%Air tank temperature 
            obj.air_on(operation)
        end
        
        function []= air_on(obj,operation)
            obj.vdot_FCV002 = operation.vdot_air_main;
            obj.vdot_FCV001 = operation.vdot_air_cooling;
            obj.vdot_FCV005 = operation.vdot_air_pilot;
            obj.vdot_LN001 = obj.vdot_FCV001 + obj.vdot_FCV002 + obj.vdot_FCV005;
            if obj.vdot_LN001<9700 && obj.vdot_LN001>obj.PCV002_highlim
                % increase cooling air if total air flow is in between
                % operation range of both pressure reducers.
                obj.vdot_FCV001 = 9700-(obj.vdot_FCV002 + obj.vdot_FCV005)+obj.PCV002_lowlim;
                obj.vdot_PCV002 = obj.PCV002_lowlim;
                obj.vdot_PCV001 = 9700;            
            elseif obj.vdot_LN001<obj.PCV002_highlim
                obj.vdot_PCV002 = obj.vdot_LN001;
                obj.vdot_PCV001 = 0;            
            else
                vdot1 = 9700;
                vdot2 = obj.vdot_LN001-vdot1;
                if vdot2 > obj.PCV002_highlim
                    vdot2 = obj.PCV002_highlim;
                    vdot1 = obj.vdot_LN001 - vdot2;
                end
                obj.vdot_PCV002 = vdot2;
                obj.vdot_PCV001 = vdot1;
            end
                
        end
        
    end
end
