


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
syms a b t
	
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

t = 0:0.1:10;
% vg = sin(2*pi*t);
% dif_vg = 2*pi*cos(2*pi*t);

vg = 0.5*t.^2;
dif_vg = t;


muy_0 = 0.1;
muy_1 = 0.5;
m = 1;
f = -muy_0*m*sign(dif_vg) - muy_1*vg;
plot(t,f)
hold on
plot(t,vg)
hold on 
plot(t, m.*dif_vg)
