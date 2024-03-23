classdef Rob < handle
    %ROB Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant = true)
      p_origin = [-0.35 0.65];
      ro_xy = Rob.robot_shape([0 0]);
      m=2;
      r = 0.08;
      K =400;
   end
    
    properties
        %local properties
        p = Rob.p_origin;
        v = [0 0 ];
        omega = 0;
        head = -pi/2;%heading
        
        %interact (wwith obs) properties
        theta = [];
        f     = [];
        r_xy = Rob.robot_shape(Rob.p_origin);
        r_heading = [0 0; 0 0];
        side;
         
    end
    
    methods (Static)
        function r_xy = robot_shape(p_origin)%position of rob
                     % khối xe
                     r = 0.08; 
                     t_s = 0:pi/180:2*pi;
                     r_x = r*cos(t_s)+p_origin(1);
                     r_y  = r*sin(t_s)+p_origin(2);
                     r_xy = [r_x;r_y];
        end
        
    end
    
    methods
        function updatePO(obj,dt) %newValue : 1x2
            %obj.p = obj.p+ obj.v*dt;
            obj.head = obj.head + obj.omega * dt;
            x = obj.p(1) + norm(obj.v)*cos(obj.head)*dt; 
            y = obj.p(2) + norm(obj.v)*sin(obj.head)*dt; 
            obj.p = [x, y];
            obj.r_xy = obj.ro_xy + obj.p';
            obj.r_heading = [obj.p; obj.p + [obj.r*cos(obj.head), obj.r*sin(obj.head)] ];
            
        end
            
        
        function plot_rob(obj)
            plot(obj.r_xy(1,:), obj.r_xy(2,:))
            hold on
            plot(obj.r_heading(:,1), obj.r_heading(:,2))
        end
        end
end
