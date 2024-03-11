


% clear all
% close all
% obs = Obs();
% a = [1 0];
% b = [1 1];


%robs = [Rob(), Rob() ];
% robs(1,2) = Rob();
% robs(2).v = [0 2];
% robs(1).p = [0 0 ];
% robs(1).p = [0 0];
% i = 0;
% v =[];
% syms x
% vo  = -1;
% %vo = -1/100100;
% f=0;
% dt = 0.1;
% m=1
% while(i<1)
%     
%     a = double(solve(m*x -f + 1*m*(vo+x*0.1)/(abs(vo+x*0.1)) +1*(vo+x*0.1) == 0,x))
%     
%     %no = double(solve(k*x^2-(k*vo+5+n)*x==0,x))
%     no = vo + a*0.1
%     no1 = (m/dt*vo -1*m)/(m/dt+1)
%     no2 = (m/dt*vo +1*m)/(m/dt+1)
%     
%     
%     warning('SOS')
%     %a = no;
%     
%     %vo = vo+1/2*double(a)*0.1^2;
%     
%     %v =[v,double(vo)];
%     i = i+1;
% end

%plot(v)
% p_x = [0 0.8 0.4 0]*2-0.4;
% p_y = [0 0 0.4 0.4]*2-0.4;
% [l_x, l_y] = make_shape(p_x,p_y);
%     min_dis = inf;
% for i = 1:length(l_x)    
%     dis = norm([0 -0.2 ]- [l_x(i) l_y(i)]);
% 
%     if dis <= min_dis 
%         min_dis = dis;
%     end
% end 
% plot(l_x, l_y)
% hold on
% axis([-5,5,-5,5])
% head = pi/3;
% pos = [3 ; 3]
%       matrix_t = [cos(head) -sin(head) pos(1) ;...
%                   sin(head) cos(head) pos(2);...
%                   0            0        1   ];
%      for i = 1:length(l_x)
%         no = (matrix_t)*[l_x(i); l_y(i);1];
%         l_x(i) = no(1);
%         l_y(i) = no(2);
%      end
%  plot(l_x, l_y)
% p1 = [0 0 ];
% p2 = [4 0 ];
% v = [2 0];

% %base = b1
% vec_dr = [p2(1) - p1(1), p2(2) - p1(2)];
% % dr_w = rotate vec_dr 90;
% matrix_t = [0 1 p1(1) ; -1 0 p1(2) ; 0 0 1];
% vec_w = inv(matrix_t)*[p2(1); p2(2); 1];
% vec_w = vec_w(1:2)' ;
% plot([0 vec_dr(1)], [0 vec_dr(2)]);
% hold on;
% plot([0 vec_w(1)], [0 vec_w(2)]);
% v_dr = dot(v, vec_dr/norm(vec_dr))*(vec_dr/norm(vec_dr));
% v_w = dot(v, vec_w/(norm(vec_w )))*(vec_w/norm(vec_w));
% if dot(v_dr, vec_dr) < 0 
%     v_dr = [0 ,0];
% end



% 
% function [m_x, m_y]=make_line(po,vec)
% 
% 
% % creat point square
% delta_l = 0.001;
% angel = atan2(vec(2), vec(1));
% delta_x = delta_l*cos(angel);
% delta_y = delta_l*sin(angel);
% 
% 
% m_x = [];
% m_y = [];
% % chia doan thang thanh nhung diem nho
% for i = 1:floor(norm(vec)/delta_l)
%     m_x = [m_x, po(1)+delta_x*i];
%     m_y = [m_y, po(2)+delta_y*i];
% end
% end
% function [s_x, s_y] = make_shape(p_x, p_y)
% % Nối điểm
% p_x = [p_x, p_x(1)];
% p_y = [p_y, p_y(1)];
% 
% s_x = [];
% s_y = [];
% 
% 
% for i = 1:length(p_x)-1
%     [l_x, l_y] = make_line([p_x(i), p_y(i)],[p_x(i+1)- p_x(i), p_y(i+1) - p_y(i)]);%tao cac diem nam gia 2 diem 
%     s_x = [s_x, l_x];
%     s_y = [s_y, l_y];
% end
% 
% 
% end

	
%f = abs(a)*abs(t);
%f = sign ()
%f_FT = fourier(f)
% syms vg(t)
% ode = -0.1*sign(diff(vg,1))*5 - 0.2*vg  == 10 ;
% cond = vg(0) == 0.5;
% ySol(t) = dsolve(ode,cond)
% S = dsolve(ode,'Implicit',true)

% s = dsolve('-(Dy)*0.1*1/abs(Dy) - 0.2*y = 10','t')
% ezplot(s,[-2,4])


% vg = sin(2*pi*t);
% dif_vg = 2*pi*cos(2*pi*t);
% Định nghĩa phương trình phi phân

% ode = @(t, v) 6 -0.25*0.67* v/norm(v)-4*v;
% 
% % Đặt điều kiện ban đầu
% initial_condition = 0.001;
% 
% % Đặt khoảng thời gian
% t_span = [0 0.1];
% 
% % Giải phương trình phi phân sử dụng ode45
% % tic
% % [t, y] = ode45(ode, t_span, initial_condition);
% % toc
% % y2 = [0];
% % for i =1:length(t)-1
% %     y2(i+1) = y2(i) + (t(i+1)-t(i))*y(i);
% % end
% % % Vẽ đồ thị kết quả
% % figure
% % plot(t, y, 'b', 'LineWidth', 1.5);
% % xlabel('Time');
% % ylabel('v(t)');
% % title('Solution of dv/dt = ...');
% % grid on;
% % figure
% % plot(t, y2, 'b', 'LineWidth', 1.5);
% % xlabel('Time');
% % ylabel('s(t)');
% % title('S');
% % grid on;
% % 
% 
% % Simulation parameters
% mass = 0.67;  % mass of the object
% friction_coefficient = 0.5;  % coefficient of friction
% total_time = 40.0;  % total simulation time
% time_step = 0.1;  % time step for simulation
% 
% % Run the simulation
% figure;
% 
% [time_values, position] = simulate_motion_with_friction_2d(mass, friction_coefficient, @applied_force, total_time, time_step);
% title('Motion with Friction and Force (2D)');
% xlabel('X Position');
% ylabel('Y Position');
% % Plot the results
% figure;
% plot(position(1, :), position(2, :));
% title('Motion with Friction and Force (2D)');
% xlabel('X Position');
% ylabel('Y Position');
% grid on;
% function [time_values, position] = simulate_motion_with_friction_2d(mass, friction_coefficient, force_function, total_time, time_step)
%     % Initialize variables
%     time_values = 0:time_step:total_time;
%     position = zeros(2, length(time_values));
%     velocity = zeros(2, length(time_values));
% 
%     % Initial conditions
%     position(:, 1) = [0; 0];
%     velocity(:, 1) = [1; 1];
%     
%     % Simulation loop
%     for i = 2:length(time_values)
%         % Calculate acceleration considering friction and applied force
%         acceleration = (force_function(time_values(i)) - friction_coefficient * velocity(:, i - 1)) / mass;
% 
%         % Update velocity and position using Euler's method
%         velocity(:, i) = velocity(:, i - 1) + acceleration * time_step;
%         position(:, i) = position(:, i - 1) + velocity(:, i) * time_step;
%         plot(position(1, :), position(2, :));
%         pause(0.01);
%         hold off;
%     end
% end
% 
% % Example force function (you can replace this with your own)
% function force = applied_force(t)
%     % A simple example of a force that changes over time
%     force = [exp(t); cos(t)]; % 2D force vector
% end



% Simulation parameters
% mass = 1.0;  % mass of the object
% friction_coefficient = 0.1;  % coefficient of friction
% total_time = 100.0;  % total simulation time
% time_step = 0.1;  % time step for simulation

% Run the simulation
% [time_values, position] = simulate_motion_with_friction_force_2d(mass, friction_coefficient, @applied_force, total_time, time_step);
% 
% % Plot the results
% figure;
% plot(position(1, :), position(2, :));
% title('Motion with Friction and Force (2D)');
% xlabel('X Position');
% ylabel('Y Position');
% grid on;
% function [time_values, position] = simulate_motion_with_friction_force_2d(mass, friction_coefficient, applied_force_function, total_time, time_step)
%     % Initialize variables
%     time_values = 0:time_step:total_time;
%     position = zeros(2, length(time_values));
%     velocity = zeros(2, length(time_values));
% 
%     % Initial conditions
%     position(:, 1) = [0; 0];
%     velocity(:, 1) = [1; 1];
% 
%     % Simulation loop using Runge-Kutta method
%     for i = 2:length(time_values)
%         % Runge-Kutta integration for each dimension
%         k1 = velocity(:, i - 1);
%         l1 = acceleration(position(:, i - 1), velocity(:, i - 1), applied_force_function(time_values(i)), mass, friction_coefficient);
% 
%         k2 = velocity(:, i - 1) + 0.5 * time_step * l1;
%         l2 = acceleration(position(:, i - 1) + 0.5 * time_step * k1, k2, applied_force_function(time_values(i)), mass, friction_coefficient);
% 
%         k3 = velocity(:, i - 1) + 0.5 * time_step * l2;
%         l3 = acceleration(position(:, i - 1) + 0.5 * time_step * k2, k3, applied_force_function(time_values(i)), mass, friction_coefficient);
% 
%         k4 = velocity(:, i - 1) + time_step * l3;
%         l4 = acceleration(position(:, i - 1) + time_step * k3, k4, applied_force_function(time_values(i)), mass, friction_coefficient);
% 
%         % Update position and velocity
%         position(:, i) = position(:, i - 1) + (1/6) * time_step * (k1 + 2*k2 + 2*k3 + k4);
%         velocity(:, i) = velocity(:, i - 1) + (1/6) * time_step * (l1 + 2*l2 + 2*l3 + l4);
%     end
% end
% 
% % Example force function (you can replace this with your own)
% function force = applied_force(t)
%     % A simple example of a force that changes over time
%     force = [sin(t)*0; cos(t)*0]; % 2D force vector
% end
% 
% % Function to calculate acceleration
% function a = acceleration(position, velocity, force, mass, friction_coefficient)
%     a = (force - friction_coefficient * velocity) / mass;
% end
close all

obs = Obs(1);
% obs.p_co = obs.so_xy(:,2002)';
obs.heading = pi/2;
% f = 5;
% obs.f = f * [cos(obs.heading), sin(obs.heading)];

robs(1,1) = Rob();
robs(1).p = [-0.4 -1];
robs(1).v = [0 0.2];
robs(1).head = pi/2;

dt = 0.1;

% get(gcf);
set(gcf,'Position',[2200 150 800 600]);

 
for i=1:200
    %obs.f = f * [cos(obs.heading), sin(obs.heading)];
    %obs.p_co = obs.s_xy(:,9100)';
    %sensing

    
    
    %apply_force if in pushing zone
    apply_Force(robs, obs, 2);
    
    robs(1).updatePO(dt)
    obs.update_aV2(dt);
    obs.updatePO(dt);
    
    

%     plot (obs.p_co(1),obs.p_co(2),'.');
%     hold on

    robs(1).plot_rob()
    

    hold on
    plot(obs.s_xy(1,:),obs.s_xy(2,:));


    axis([-5,5,-5,5])
    grid on


    pause(0.1);
    hold off

    
end


