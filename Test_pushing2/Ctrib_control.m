classdef Ctrib_control <handle
    %CONTROL Summary of this class goes here
    % control outer loop (contributed) , cal f_fid for each robot with
    % input is error angles beetween designed direct and v robot (alpha)
    
    %   Detailed explanation goes here
    % step 1: set up side for each robot (not automatic for this prj) ->
    % know nL, nR
    % step 2: apply PID to cal f_id
    
    properties
        
        f0 = 1.5;
        K_Po = 2;
        K_Io = 0.5;
        K_Do = 0;
        error_old = []; %(alpha) dau vao dieu khien tai thoi diem (t-dt)
        f_di  = []; % dau ra dieu khien
        nL=0;
        nR=0;
        f_max = 4;
        
    end
    
    methods
        
        function obj = Ctrib_control(robs)
            if nargin ==1
                N = length(robs);
                obj.error_old = zeros(1,N);
                obj.f_di = zeros(1,N);
               for i = 1:N
                   if robs(i).side == 1
                       obj.nR = obj.nR+1;
                   else
                       obj.nL = obj.nL+1;
                   end
                end

            end

        end
        
        function apply(obj, robs, designed_angle,dt)
            num = [obj.nL 0 obj.nR];
            for i=1:length(robs)
                %sensing and cal error (alpha at time t)
                alpha = designed_angle - (robs(i).head - floor((robs(i).head + pi)/(2*pi))*2*pi ); %%new_eror
                %chuan hoa head thuoc [-pi,pi];
                propor = obj.K_Po * alpha;
                inter  = obj.error_old(i) + alpha*dt;
                deriva = obj.K_Do*alpha/dt;
                
                side = robs(i).side;
                f_dii = (obj.f0 + side*(propor + inter + deriva))/num(side+2);
                %chuan hoa
                if f_dii < 0 
                    obj.f_di(i) = 0;
                elseif f_dii > obj.f_max
                    obj.f_di(i) = obj.f_max;
                else
                    obj.f_di(i) = f_dii;
                end
                
            end
        end
    end
end

