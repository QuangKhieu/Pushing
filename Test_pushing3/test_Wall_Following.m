close all

obs = Obs(1);
% obs.p_co = obs.so_xy(:,2002)';
obs.heading = pi/2;

% f = 5;
% obs.f = f * [cos(obs.heading), sin(obs.heading)];

robs(1,4) = Rob();

robs(1).p = [-0.5 -1];
robs(1).v = [0 0.1];
robs(1).head = pi/2;

robs(2).p = [-0.2 -1];
robs(2).v = [0 0.1];
robs(2).head = pi/2;

robs(3).p = [0.2 -1];
robs(3).v = [0 0.1];
robs(3).head = pi/2;

robs(4).p = [0.5 -1];
robs(4).v = [0 0.1];
robs(4).head = pi/2;

robs(1).p = [-0.4 -1];
robs(1).v = [0 0.1];
robs(1).head = pi/2;





dt = 0.1;

% get(gcf);
set(gcf,'Position',[2200 150 800 600]);


for i=1:600
    %obs.f = f * [cos(obs.heading), sin(obs.heading)];
    %obs.p_co = obs.s_xy(:,9100)';
    %sensing


   w_f(robs, obs)

    
    %apply_force if in pushing zone
    %apply_Force(robs, obs, 0);
    
    
    %obs.update_aV2(dt);
    %obs.updatePO(dt);
    
    

     %plot (obs.p_co(1),obs.p_co(2),'.');
        for k =1:length(robs)
            robs(k).updatePO2(dt)
            robs(k).plot_rob()
            hold on
        end


    
    

    

    plot(obs.s_xy(1,:),obs.s_xy(2,:));

    axis([-2,2,-2,2])
    grid on

    pause(0.1);
    hold off

    
end
function  w_f(robs, obs)
%ref_obj :(2xN), obj(1x2) point colision
    n_clear = [0, 1];   
    error = 0;
    rS = 0.12;
    e1 = 0;
    e2 = 0.1 ;
    e3 = 0.9;
    e4 = 1.3;
    
    p1 = 0.1;
    p2 = 48;
    
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
        %dis = min_dis
            if min_dis <= 0.2
                check = 1;
            else
                check = 0;
            end
            
            %check_zone = check;

            obs.p_co(i_th,:) = p_co;
            %pc = p_co;
            if (check == 1) %neu va cham tinh toan huong , do lon cua luc tac dung vao objs
                delta_i_ = p_co - robs(i_th).p;
                delta_i = delta_i_/norm(delta_i_)*(0.08-norm(delta_i_));
                f_ = robs(i_th).K*delta_i;
               % obs.f(i_th,:) = f_;

                robs(i_th).theta = atan2(delta_i(2), delta_i(1))...
                                    -(robs(i_th).head - floor((robs(i_th).head + pi)/(2*pi))*2*pi );
                %robot nhan luc nguoc chieu, cung do lon
                robs(i_th).f  = -1*f_;   
                %a = robs(i_th).theta;
                
                n_x = norm(delta_i_); % norm of x_ti
                % cal q_, phi
                q_ = (n_x/rS - e3)/e4;
                
                
                phi = e1 + e2* (q_/1 + abs(q_) );
                
                v_odm = phi * (delta_i_)/n_x;
                
                % cac truong hop
                if i_th ~=4
                 
                    x_ = [-delta_i_(2), delta_i_(1) ];
                    
                else
                    x_ = [delta_i_(2), delta_i_(1) ];
                end
                
                 error = rad2deg(atan2(n_clear(2), n_clear(1)) ...
                -   atan2(delta_i_(2), delta_i_(1)));
                


                v_bf = (e1+e2/2)*x_/norm(x_);
                if i_th == 4 && error < -90
                       
                        v_bf = 0;
                end
                if i_th == 1 && error > 90
                        v_bf = 0;
                end
%                 if i_th ==4
%                   error
%                 end
                %distance
                d1 = norm (robs(1).p - robs(2).p);
                d2 = norm (robs(2).p - robs(3).p);
                d3 = norm (robs(3).p - robs(4).p);
                if i_th == 2
                n_x = norm(x_); % norm of x_ti
                % cal q_, phi

                    q_ = (d1 - d2)/(d1+d2)


                    phi = -p1/2 + p1/(1 + exp(-(p2*q_)))

                    v_bf = phi * (x_)/n_x     ;            
                end
                
                if i_th == 3
                    q_ = (d2 - d3)/(d2+d3);


                    phi = -p1/2 + p1/(1 + exp(-(p2*q_)));

                    v_bf = phi * (x_)/n_x;                     
                end
                
                robs(i_th).v = v_odm + v_bf;


                
            else%neu khong va cham
               %obs.f(i_th,:) = [0 0]; 
               robs(i_th).theta = inf;
               robs(i_th).f     = [0 0];
            end
            
            %%% set flag to warning rob too near obs or too far from obs
            % too far: continue inner control (PI controler) using v w of
            % previous looptime(t-deltat)
            % too near: v,w = 0 
            
    end

end

