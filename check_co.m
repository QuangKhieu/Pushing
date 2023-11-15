function check_co(robs, obs)
%ref_obj :(2xN), obj(1x2) point colision
       
    ref_obj = obs.s_xy;
    % xet tung robot
    for i_th = 1:length(robs)
        min_dis = inf;
        obj = robs(i_th).p;
        for i = 1:length(ref_obj)    
            dis = norm(obj-ref_obj(:,i)');
            if dis <= min_dis 
                min_dis = dis;
                p_co = ref_obj(:,i)';
            end
        end  

            if min_dis <= 0.08
                check = 1;
            else
                check = 0;
            end

            obs.p_co(i_th,:) = p_co;
            if (check == 1)
                delta_i_ = p_co - robs(i_th).p;
                delta_i = delta_i_/norm(delta_i_)*(0.08-norm(delta_i_));
                obs.f(i_th,:) = robs(i_th).K*delta_i;
            else
               obs.f(i_th,:) = [0 0]; 
            end
    end
end