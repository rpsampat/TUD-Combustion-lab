classdef Operation<handle
    properties
        vdot_air_main%lnpm
        vdot_air_cooling%lnpm
        vdot_air_pilot%lnpm
    end
    methods
        function obj = Operation(settings, main_burner, comb, exhaust)            
            mdot_gas = main_burner.mdot_air + main_burner.mdot_fuel;        
            obj.vdot_air_main = main_burner.vdot_air;
            obj.vdot_air_cooling = exhaust.vdot_air_cooling;
            obj.vdot_air_pilot = 0;            
        end
    end
end
