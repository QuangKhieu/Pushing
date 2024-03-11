function apply_Force(robs, obs, f)
%APPLY_FORCE Summary of this function goes here
% apply wanted_force to obs when robot go into pushing zone  
%   Detailed explanation goes here
%ref_obj :(2xN), obj(1x2) point colision
       
%    ref_obj = obs.s_xy;
    
%     % xet tung robot
     for i_th = 1:length(robs)
%         min_dis = inf;
%         obj = robs(i_th).p;
%         for i = 1:length(ref_obj)    
%             dis = norm(obj-ref_obj(:,i)');
%             if dis <= min_dis 
%                 min_dis = dis;
%                 p_co = ref_obj(:,i)';
%             end
%         end  
% 
%             if min_dis <= 0.08
%                 check = 1;
%             else
%                 check = 0;
%             end
% 
%             obs.p_co(i_th,:) = p_co;
%             if check == 1
%                 obs.f = f * [cos(obs.heading), sin(obs.heading)];               
%             else
%                 obs.f = [0, 0];
%             end
        if i_th ~=1 || i_th~=4
            obs.f(i_th, :) = f * (robs(i_th).f)/norm(robs(i_th).f);
        else
             obs.f(i_th, :) = [0 0];
        end
    end
end

