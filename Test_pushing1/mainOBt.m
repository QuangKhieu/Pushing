close all
clear classes
count = 0;
obs = Obs();
rob = Rob();
dt =0.1;
rob.v = [0 -0.4];
obs.p = [0 0.7];
vObs=[];


while (1)
    rob.v = [0 -0.3];
    [check,p_co] = check_co(rob,obs);
    if check == 1 && (norm(rob.v) > norm(obs.v))
        count = count + 1;
        v2 = obs.v - 2*rob.m/(rob.m+obs.m)*(dot(obs.v - rob.v, p_co - rob.p))...
                /(norm(p_co - rob.p))^2*(p_co - rob.p);
        v1 = rob.v - 2*obs.m/(rob.m+obs.m)*(dot(rob.v - obs.v, rob.p - p_co))...
                /(norm(rob.p - p_co))^2*(rob.p - p_co);  
        %obs.v = v2;
        obs.updateV(v2,dt);

         rob.v = dotVec(obs.v,rob.v);
                 rob.v = v1;
%                  rob.v = [0 0];
        obs.vec_ana(p_co); 
%         if count ==1
%         break
%         end
    end

    obs.updateV([],dt);
    if check && (norm(rob.v) > norm(obs.v))
       % rob.v = dotVec(obs.v,[0.2 -0.6]);
    end
    rob.updatePO(dt);
    obs.updatePO(dt);
    
    if ~check
    vObs = [vObs, norm(obs.v)];
    hold off
    plot(obs.s_xy(1,:),obs.s_xy(2,:))
    hold on
    plot(rob.r_xy(1,:),rob.r_xy(2,:))
    axis([-5,5,-5,5])
    pause(0.01)
    end

%         if count ==1
%         break
%         end   
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
function c = dotVec(a,b) %chiếu a lên b
    c = dot(a,b/norm(b))*b/norm(b);
end