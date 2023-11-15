clear

p_x = [0 0.4 0.4 0]*2-0.4;
p_y = [0 0 0.4 0.4]*2-0.4;

[l_x, l_y] = make_shape(p_x,p_y);
l_xn = l_x;
l_yn = l_y;
[r_x, r_y] = robot_shape();


obs = Obs();
rob = Rob();
rob.v = [0 -0.4];
obs.p = [-4 -0.2];
dt =0.1;
r_x = r_x + rob.p(1);
r_y = r_y + rob.p(2);
p_co = obs.p;


%%
count = 0;
pobs=[0];
v_dr = [0 0];
v_w = [0 0];
muy = 0; %he so ma sat
fms = muy*obs.m*9.8;
v_ms = fms*dt/obs.m;
obs.v_ms = v_ms;
while(1)
   % rob.v = [0.4 0];
    
    %[check, p_co] = check_colision([rob.p(1) ,rob.p(2) ],[l_xn', l_yn']);%p_co is point_colision
    [check,p_co] = check_co(rob,obs);
    %p_co = obs.p;
    if check == 1 && (norm(rob.v) > norm(obs.v))
%
        count = count + 1;
%         no = inv([rob.m obs.m; rob.m -obs.m])*[rob.m*rob.v(1)+obs.m*(obs.v(1));0];
%         obs.v(1) = no(2);
%         no = inv([rob.m obs.m; rob.m -obs.m])*[rob.m*rob.v(2)+obs.m*(obs.v(2));0];
%         obs.v(2) = no(2);
        
        v2 = obs.v - 2*rob.m/(rob.m+obs.m)*(dot(obs.v - rob.v, p_co - rob.p))...
                /(norm(p_co - rob.p))^2*(p_co - rob.p);
        v1 = rob.v - 2*obs.m/(rob.m+obs.m)*(dot(rob.v - obs.v, rob.p - p_co))...
                /(norm(rob.p - p_co))^2*(rob.p - p_co);  
%         v2 = obs.v - 2*rob.m/(rob.m+obs.m)*(dot(obs.v - rob.v, obs.p - rob.p))...
%                 /(norm(obs.p - rob.p))^2*(obs.p - rob.p);
%         v1 = rob.v - 2*obs.m/(rob.m+obs.m)*(dot(rob.v - obs.v, rob.p - obs.p))...
%                 /(norm(rob.p - obs.p))^2*(rob.p - obs.p);
        %obs.v = v1;
        

        obs.v = v2 - [-v_ms 0];
        rob.v = v1;
        obs.vec_ana(p_co);
% 
%         break
%         [check,~,v_dr, v_w] = vec_ana(p_co, obs.p , obs.v, obs.v_ms);
%         
%         r = norm(obs.p - p_co);
%         obs.omega = check*norm(v_w)/r;
%         break
        %rob.v =obs.v;%???
    end
    
%     if ~check && norm(obs.v)~=0
%         [check,obs.v,v_dr, v_w] = vec_ana(p_co, obs.p , obs.v, obs.v_ms);
%         obs.omega = check*norm(v_w)/r;
%     end
    pobs=[pobs,(norm(obs.v))];

%     
%     if(obs.v(2)>0.25)
%         obs.v(2) = obs.v(2) - 0.25;
%         
%     end

%     obs.heading = obs.heading + obs.omega*dt;
%     obs.p = obs.p +[obs.v_dr(1)*dt, obs.v_dr(2)*dt];
    rob.updatePO(dt);
    obs.updatePO(dt);
    hold off
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    plot(rob.r_xy(1,:),rob.r_xy(2,:))
    %[v_dr, v_w] = vec_ana(p_co, obs.p , obs.v, obs.v_ms);
    
%     l_x = l_x + obs.v(1)*dt;
%     l_y = l_y + obs.v(2)*dt;

%         rob.p = rob.p +[rob.v(1)*dt, rob.v(2)*dt];
%         r_x  = r_x + rob.v(1)*dt;
%         r_y  = r_x + rob.v(2)*dt;
% 
% 
%     plot(r_x, r_y);
    axis([-5,5,-5,5])
    pause(0.01)
   
   % [l_xn, l_yn]= draw(l_x, l_y, r_x, r_y, obs.heading, obs.p);  
    
    
end




%% 
for i =0:0
   
     hold off
     plot(l_x, l_y)
     hold on
    check = check_colision([0.4 ,-1+0.1*i ],[l_x', l_y']);
    plot(r_x+0.4, r_y-1+0.1*i)
    axis([-3,3,-3,3])
    pause(0.10)
   
end
%%
function [check,v_update, v_dr, v_w] = vec_ana(p1, p2 , v, v_ms)%p1 la diem va cham, p2 la obs.p

if norm(v) <= v_ms
    v = [0 0 ];
else
v_ms = [ v_ms*cos(atan2(-v(2), -v(1))) ,v_ms*sin(atan2(-v(2), -v(1))) ];
v = v + v_ms;
end
v_update = v;
%% Chiếu vector có đầu vào gồm 2 điểm 1 vector V (p1 p2 vc)
%base = p1
vec_dr = [p2(1) - p1(1), p2(2) - p1(2)];
% dr_w = rotate vec_dr 90;
matrix_t = [0 1 p1(1) ; -1 0 p1(2) ; 0 0 1];
vec_w = inv(matrix_t)*[p2(1); p2(2); 1];
vec_w = vec_w(1:2)' ;

% phi =pi/2;
% x_ = p1(1) + (p2(1) - p1(1))*cos(pi/2) - (p2(2) - p1(2))*sin(phi);
% y_ = p1(2) + (p2(1) - p1(1))*sin(pi/2) + (p2(2) + p1(2))*cos(phi);

%Chiếu vector
v_dr = dot(v, vec_dr/norm(vec_dr))*(vec_dr/norm(vec_dr));
v_w = dot(v, vec_w/(norm(vec_w )))*(vec_w/norm(vec_w));

cro = cross([p1 - p2,0],[v_w,0]);
if cro(3)>0 %check direct of omega
    check = 1;
else
    check = -1;
end
if dot(v_dr, vec_dr) < 0 
    v_dr = [0 ,0];
end


end
function [l_xn, l_yn] = draw(l_x, l_y, r_x, r_y, head, pos)

      hold off
      matrix_t = [cos(head) -sin(head) pos(1) ;...
                  sin(head) cos(head) pos(2);...
                  0            0        1   ];
     for i = 1:length(l_x)
        no = (matrix_t)*[l_x(i); l_y(i);1];
        l_xn(i) = no(1);
        l_yn(i) = no(2);
     end
     
    plot(l_xn, l_yn)
    hold on
    plot(r_x, r_y)

    axis([-5,5,-5,5])
    pause(0.01)
    
end

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
function [s_x, s_y] = make_shape(p_x, p_y)
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


end
function [check, p_co] = check_colision(obj, ref_obj) %ref_obj :(Nx2), obj(1x2) point colision
    min_dis = inf;
for i = 1:length(ref_obj)    
    dis = norm(obj-ref_obj(i,:));

    if dis <= min_dis 
        min_dis = dis;
        p_co = ref_obj(i,:);
    end
end  


    if min_dis <= 0.13
        check = 1;
    else
        check = 0;
    end

end
function [r_x, r_y]=robot_shape()
% khối xe
 r = 0.1; 
 t_s = 0:pi/180:2*pi;
 r_x = r*cos(t_s);
 r_y  = r*sin(t_s);
end
function [check, p_co] = check_co(rob, obs)
%ref_obj :(2xN), obj(1x2) point colision
obj = rob.p;
ref_obj = obs.s_xy;
    min_dis = inf;
for i = 1:length(ref_obj)    
    dis = norm(obj-ref_obj(:,i)');
    if dis <= min_dis 
        min_dis = dis;
        p_co = ref_obj(:,i)';
    end
end  


    if min_dis <= 0.13
        check = 1;
    else
        check = 0;
    end





end