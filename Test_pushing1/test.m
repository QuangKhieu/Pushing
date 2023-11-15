
obs =Obs();
rob  =Rob();
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
    disp(min_dis)

    if min_dis <= 0.13
        check = 1;
    else
        check = 0;
    end

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
%% Chiếu vector có đầu vào gồm 2 điểm 1 vector V (p1 p2 vc)
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