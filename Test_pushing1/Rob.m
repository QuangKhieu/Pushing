classdef Rob < handle
    %ROB Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant = true)
      p_origin = [0 2];
      ro_xy = Rob.robot_shape([0 0]);
      m=2;
   end
    
    properties

        p = Rob.p_origin;
        v = [0 0 ];
        r_xy = Rob.robot_shape(Rob.p_origin);

         
    end
    
    methods (Static)
        function r_xy = robot_shape(p_origin)%position of rob
                     % khối xe
                     r = 0.1; 
                     t_s = 0:pi/180:2*pi;
                     r_x = r*cos(t_s)+p_origin(1);
                     r_y  = r*sin(t_s)+p_origin(2);
                     r_xy = [r_x;r_y];
        end
        
    end
    
    methods
        function updatePO(obj,dt) %newValue : 1x2
            obj.p = obj.p+ obj.v*dt;
            obj.r_xy = obj.ro_xy + obj.p'; 
        end
            
        

    end
end
