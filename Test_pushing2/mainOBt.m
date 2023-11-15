close all
clear classes
count = 0;

robs(1,2) = Rob();
obs = Obs(2);
robs(2).p = [0.35 0.65];
robs(1).v = [0 -0.05];
robs(2).v = [0 -0.04];
dt =0.1;
k = 50;
figure(1)
trace = [obs.p];
v = [obs.v];



while(true)
    
     check_co(robs,obs);

%     obs.f = k*[(obs.s_xy(1,1700) - obs.s_xy(1,1600)) , (obs.s_xy(2,1700) - obs.s_xy(2,1600))];
%     %obs.p_co = [obs.s_xy(1,166), obs.s_xy(2,166)];
%    obs.p_co = [obs.s_xy(1,1333), obs.s_xy(2,1333)];


    obs.update_aV(dt);
    obs.updatePO(dt);
    robs(1).updatePO(dt);
    robs(2).updatePO(dt);
    v = [v, norm(obs.v)];
    
    hold off
   
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    plot(robs(1).r_xy(1,:), robs(1).r_xy(2,:))
    hold on
    plot(robs(2).r_xy(1,:), robs(2).r_xy(2,:))
    
    trace = [trace; obs.p];
    plot(trace(:,1),trace(:,2))
    axis([-3,3,-3,3])
      count = count +1;
      
    pause(0.01);
     
end

