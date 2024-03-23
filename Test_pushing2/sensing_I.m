function sensing_I(robs, obs)
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
            if (check == 1) %neu va cham tinh toan huong , do lon cua luc tac dung vao objs
                delta_i_ = p_co - robs(i_th).p;
                delta_i = delta_i_/norm(delta_i_)*(0.08-norm(delta_i_));
                f_ = robs(i_th).K*delta_i;
                obs.f(i_th,:) = f_;

                robs(i_th).theta = atan2(delta_i(2), delta_i(1))...
                                    -(robs(i_th).head - floor((robs(i_th).head + pi)/(2*pi))*2*pi );
                %robot nhan luc nguoc chieu, cung do lon
                robs(i_th).f  = -1*f_;   
                a = robs(i_th).theta
            else%neu khong va cham
               obs.f(i_th,:) = [0 0]; 
               robs(i_th).theta = inf;
               robs(i_th).f     = [0 0];
            end
            
            %%% set flag to warning rob too near obs or too far from obs
            % too far: continue inner control (PI controler) using v w of
            % previous looptime(t-deltat)
            % too near: v,w = 0 
            
    end
end