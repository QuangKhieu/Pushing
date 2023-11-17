close all
clear 
count = 0;
%Setup 2 robot, 1 object
robs(1,2) = Rob();
obs = Obs(2);
robs(1).p = [-0.25 0.55];
robs(2).p = [0.25 0.55];
robs(1).v = [0 -0.02];
robs(2).v = [0 -0.02];
robs(1).side =  1; 
robs(2).side =  -1;
dt =0.1;
cb_control = Ctrib_control(robs);
in_control =  Inner_control(robs);
k = 50;
figure(1)
trace = [obs.p];
v = [obs.v];
% các ràng buộc : vật luôn tương tác với vật, góc 

designed_angle = -pi/2 + pi/18;


while(true)
    %sensing outer loop  and
    %contributed control(control outer)....
    cb_control.apply(robs, designed_angle, dt);
    
    % sensing inner loops
    sensing_I(robs,obs); 
    
    %control inner loops.... 
    f_di = cb_control.f_di;
    in_control.apply(robs, f_di, dt);
    
    
    % operate inner
    robs(1).updatePO(dt);
    robs(2).updatePO(dt);
    % operate outer
    sensing_I(robs,obs); 
    obs.update_aV(dt);
    obs.updatePO(dt);





    %v = [v, norm(obs.f)];
    
    hold off
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    robs(1).plot_rob()
    hold on
    robs(2).plot_rob()

    
    trace = [trace; obs.p];
    plot(trace(:,1),trace(:,2))
    axis([-3,3,-3,3])
      count = count +1;
      
    pause(0.01);
     
end

