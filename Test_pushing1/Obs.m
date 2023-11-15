classdef Obs < handle
    %OBS Summary of this class goes here
    %   Detailed explanation goes here
   properties (Constant = true)
       %dỉnh đa giác
        p_x = [0 0.4 0.4 0]*2-0.4; 
        p_y = [0 0 0.4 0.4]*2-0.4;
        so_xy = Obs.make_shape(Obs.p_x, Obs.p_y); % 2*N
        env = Env() ;
        m = 5;
   end
    properties

        p=[0 0];
        v= [0 0];  
        heading = 0;
        omega = 0;
        v_ms = 0;
        v_dr = [0 0];
        v_w = [0 0];
        s_xy = Obs.make_shape(Obs.p_x, Obs.p_y); % 2*N
        fms = Obs.env.muy * Obs.env.g * Obs.m;
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
                delta_l = 0.001;
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

        % phân tích lực khi va chạm Khi va chạm
        
        function vec_ana(obj,p_co)
            %p1 la diem va cham, p2 la obs.p
            p1 = p_co;
            p2 = obj.p;
            %% Chiếu vector có đầu vào gồm 2 điểm 1 vector V (p1 p2 vc)
            %base = p1
            vec_dr = [p2(1) - p1(1), p2(2) - p1(2)];
            % vec_w = rotate vec_dr 90;
            matrix_t = [0 1 p1(1) ; -1 0 p1(2) ; 0 0 1];
            vec_w = inv(matrix_t)*[p2(1); p2(2); 1];
            vec_w = vec_w(1:2)' ;

            %Chiếu vector
            obj.v_dr = dot(obj.v, vec_dr/norm(vec_dr))*(vec_dr/norm(vec_dr));
            obj.v_w = dot(obj.v, vec_w/(norm(vec_w )))*(vec_w/norm(vec_w));

            cro = cross([p1 - p2,0],[obj.v_w,0]);
            if cro(3)>0 %check direct of omega
                obj.omega = norm(obj.v_w)/norm(p2-p1);
            else
                obj.omega = -norm(obj.v_w)/norm(p2-p1);
            end
            if dot(obj.v_dr, vec_dr) < 0 
                obj.v_dr = [0 ,0];
                obj.v_w  = [0,0];
            end

        end
        %% update v
        function updateV(obj, v_up, dt)
            %khi va cham v_up :v_obs_new
            if ~isempty(v_up)
                obj.v = v_up;
                return
            end
            %khi khong va cham
            if isempty(v_up)
                v_ms = obj.fms*dt/obj.m;
                if norm(obj.v)> v_ms
                    angle = atan2(obj.v(2), obj.v(1));
                    k = v_ms/norm(obj.v);
                    obj.v(1) = obj.v(1) - v_ms*cos(angle);
                    obj.v(2) = obj.v(2) - v_ms*sin(angle);
                    obj.v_dr = obj.v_dr - k*obj.v_dr;
                    obj.v_w  = obj.v_dr - k*obj.v_w;
                else
                    obj.v = [0 0];
                    obj.v_dr = [0 0];
                    obj.v_w =[0 0];
                    obj.omega = 0;
                end
            
            end
        end
        %update p 
        function updatePO(obj, dt) 
            %cal PO of center obs
            obj.p = obj.p + obj.v_dr*dt;
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

