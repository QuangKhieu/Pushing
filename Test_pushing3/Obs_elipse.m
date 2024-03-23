classdef Obs_elipse < handle
    %OBS Summary of this class goes here
    %obs elipse shape
    %   Detailed explanation goes here
   properties (Constant = true)
        p_origin = [0 0];
        so_xy = Obs_elipse.make_shape([0, 0]); % 2*N
        env = Env() ;
        m = 2;
        I = 1;
   end
    properties
        
        p= Obs_elipse.p_origin;
        v= [0 0];  
        heading = 0;
        omega = 0;
        f = [0 0];%(m*2) m hang 2 cot
        a = [0 0];
        a_ht = 0; %d(omega)/dt
        s_xy = Obs_elipse.so_xy; % 2*N
        p_co;

    end
    
    methods (Static)
        %Khối obs
        function s_xy = make_shape(p_origin)
             a = 2;
             b = 1;
             t_s = 0:2*pi/9600:2*pi;
             s_x = a*cos(t_s)+p_origin(1);
             s_y  = b*sin(t_s)+p_origin(2);
             s_xy = [s_x;s_y];            
        end
    end
    
    methods
        function obj = Obs(number)
            if nargin == 1
                obj.f = zeros(number,2);
            end
            
        end

       
       
        function update_aV(obj, dt)
           %cal sum F 
            sum_f = [sum(obj.f(:,1)),sum(obj.f(:,2))]; % tong hop luc
 
            %update V
            if norm(obj.v) ~= 0
               % cal a
               obj.a = (sum_f - obj.env.muy0_t * obj.m*obj.v/norm(obj.v) ...
                  -  obj.env.muy1_t*obj.v) / obj.m;
               v_ = obj.v + obj.a * dt;
               
                
            else
               obj.a = sum_f / obj.m;
               v_ = obj.v + obj.a * dt;                
                
            end
            
            %cali
            
            if (v_(1)*obj.v(1) < 0)
                obj.v(1) = 0;
            else
                obj.v(1) = v_(1);
            end
            
            if (v_(2)*obj.v(2) < 0)
                obj.v(2) = 0;
            else
                obj.v(2) = v_(2);
            end
            
            %update omega
            %sum moment
            
            r_vec = obj.p_co-obj.p ;%n*2 (n la so robot)
            % cross matrix f vs matrix r
                  
            cros = cross([r_vec'; zeros(1, height(r_vec))],...
                [obj.f' ; zeros(1,height(r_vec))]);
            
            sum_cf = sum(cros(3,:))  ;  
            norm(obj.omega);
            if norm(obj.omega) ~= 0
                obj.a_ht = (sum_cf...
                    - obj.env.muy0_r*obj.I*obj.omega/norm(obj.omega)...
                    - obj.env.muy1_r*obj.omega) / obj.I;
                omega_ = obj.omega + obj.a_ht * dt;
            else
                obj.a_ht = sum_cf / obj.I;
                omega_ = obj.omega + obj.a_ht * dt;
            end
%             obj.omega = omega_;
%             
            
            if omega_*obj.omega < 0
                obj.omega = 0;
            else
                obj.omega = omega_;
            end
            
        end
        %update p 
        function updatePO(obj, dt) 
            %cal PO of center obs
            obj.p = obj.p + obj.v*dt;
            obj.heading = obj.heading + obj.omega*dt;
            
            %update PO of eadge
            head = obj.heading;
            pos = obj.p;
            matrix_t = [cos(head) -sin(head) pos(1) ;...
                          sin(head) cos(head) pos(2);...
                          0            0        1   ];
             for i = 1:length(obj.so_xy)
                no = (matrix_t)*[obj.so_xy(1,i); obj.so_xy(2,i);1];
                obj.s_xy(1,i) = no(1);
                obj.s_xy(2,i) = no(2);
             end
        end


    end
end

