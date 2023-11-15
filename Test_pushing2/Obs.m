classdef Obs < handle
    %OBS Summary of this class goes here
    %   Detailed explanation goes here
   properties (Constant = true)
       %dỉnh đa giác
        p_x = [0 0.4 0.4 0]*2-0.4; 
        p_y = [0 0 0.4 0.4]*2-0.4;
        so_xy = Obs.make_shape(Obs.p_x, Obs.p_y); % 2*N
        env = Env() ;
        m = 1;
        I = 0.05;
   end
    properties
        
        p=[0 0];
        v= [0 0];  
        heading = 0;
        omega = 0;
        f = [0 0];%(m*2) m hang 2 cot
        a = [0 0];
        a_ht = 0; %d(omega)/dt
        s_xy = Obs.make_shape(Obs.p_x, Obs.p_y); % 2*N
        p_co;

    end
    
    methods (Static)
        %Khối obs
        function s_xy = make_shape(p_x, p_y)
            % Nối điểm
            p_x = [p_x, p_x(1)];
            p_y = [p_y, p_y(1)];

            s_x = [];
            s_y = [];


            for i = 1:length(p_x)-1
                [l_x, l_y] = make_line([p_x(i), p_y(i)],[p_x(i+1)- p_x(i), p_y(i+1) - p_y(i)]);%tao cac diem nam gia 2 diem 
                s_x = [s_x, l_x];
                s_y = [s_y, l_y];
            end
            s_xy = [s_x;s_y];
            function [m_x, m_y]=make_line(po,vec)

                % creat point square
                delta_l = 0.0005;
                angel = atan2(vec(2), vec(1));
                delta_x = delta_l*cos(angel);
                delta_y = delta_l*sin(angel);


                m_x = [];
                m_y = [];
                % chia doan thang thanh nhung diem nho
                for i = 1:floor(norm(vec)/delta_l)
                    m_x = [m_x, po(1)+delta_x*i];
                    m_y = [m_y, po(2)+delta_y*i];
                end
            end
        end
    end
    
    methods
        function obj = Obs(number)
            if nargin == 1
                obj.f = zeros(number,2);
            end
            
        end

        function update_aV(obj, dt)
            % pt tong hop 
            % ma = sigma(f) -muy0_t*sign((vo + a*dt) - muy1_t*(vo + a*dt))=>solve
            % Ia_ht = sigma(cross(fi,ri)) +Mf
            
            sum_f = [sum(obj.f(:,1)),sum(obj.f(:,2))]; % tong hop luc
            
            %V_x
            v_x = 0;
            no1 = (obj.m/dt*obj.v(1) + sum_f(1) - obj.env.muy0_t*obj.m)/(obj.m/dt+obj.env.muy1_t);% v_x1
            no2 = (obj.m/dt*obj.v(1) + sum_f(1) + obj.env.muy0_t*obj.m)/(obj.m/dt+obj.env.muy1_t);% v_x2
            
            if (no1>0 && no2>0)
                v_x = no1;
            elseif (no1<0 && no2<0)
                v_x = no2;
            elseif ((no1<0 && no2>0) || no1*no2 == 0 )
                v_x = 0;
            else
                warning('problem with no at V_x')
                disp(no1,no2)
            end
            
            %V_y
            v_y = 0;
            no1 = (obj.m/dt*obj.v(2) + sum_f(2) - obj.env.muy0_t*obj.m)/(obj.m/dt+obj.env.muy1_t);% v_y1
            no2 = (obj.m/dt*obj.v(2) + sum_f(2) + obj.env.muy0_t*obj.m)/(obj.m/dt+obj.env.muy1_t);% v_Y2
            
            if (no1>0 && no2>0)
                v_y = no1;
            elseif (no1<0 && no2<0)
                v_y = no2;
            elseif ((no1<0 && no2>0) || no1*no2 == 0 )
                v_y = 0;
            else
                warning('problem with no at V_y')
                disp(no1,no2)
            end
            
            obj.a = ([v_x, v_y] - obj.v)/dt;
            obj.v = [v_x, v_y];
            
            
            
            
            
            
            
            %sum moment
            sum_cf = 0;
            r_vec = obj.p_co-obj.p; %n*2 (n la so robot)
            % cross matrix f vs matrix r
            
            
            cros = cross([r_vec'; zeros(1, height(r_vec))] , [obj.f' ; zeros(1,height(r_vec))]);
            omega_ = 0;
            sum_cf = sum(cros(3,:)); 
            
            
            no1 = (obj.I/dt*obj.omega + sum_cf - obj.env.muy0_r*obj.I)/(obj.I/dt + obj.env.muy1_r);
            no2 = (obj.I/dt*obj.omega + sum_cf + obj.env.muy0_r*obj.I)/(obj.I/dt + obj.env.muy1_r);
            if (no1>0 && no2>0)
                omega_ = no1;
            elseif (no1<0 && no2<0)
                omega_ = no2;
            elseif ((no1<0 && no2>0) || no1*no2 == 0 )
                omega_ = 0;
            else
                warning('problem with no at oemga')
                disp(no1,no2)
            end
            obj.a_ht = (obj.omega - omega_)/dt; 
            obj.omega = omega_;
            
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

