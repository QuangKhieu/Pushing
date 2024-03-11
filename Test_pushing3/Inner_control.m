classdef Inner_control < handle
    %INNER_CONTROL Summary of this class goes here
    % inner control to maining pushing, exert desinged force f_id
    %   Detailed explanation goes here
    
    properties
        K_omega = 0.4;
        K_Pin = 0.03;
        K_Iin = 0.01;
        old_error = []; %old (e_fi)
    end
    
    methods
        function obj = Inner_control(robs)
            if nargin == 1
                N = length(robs);
                obj.old_error = zeros(1,N);    
            end
        end
        
        function apply(obj, robs, f_di, dt)
            for i = 1:length(robs)
                % check f_i == 0? ==1?...
                f_i = norm(robs(i).f);
                
                if f_i == 0
                    if i==2
                        disp("vkl")
                    end
                    %robs(i).omega = 0;
                    %robs(i).omega = obj.K_omega*robs(i).theta;
                    continue;       
                end
                
                e_fi = f_di(i) - f_i ;
                
                %PI control
                
%                 a = obj.K_omega;
                b = robs(i).theta;


                robs(i).omega = obj.K_omega*robs(i).theta;
                    
                v     = obj.K_Pin*e_fi + ...
                                obj.K_Iin*(obj.old_error(i) + e_fi*dt);
                %chuan hoa  
                if v > 0.07
                     robs(i).v = 0.07;
                elseif(v<0)
                    robs(i).v = 0.02;
                else
                    robs(i).v = v;
                end
                    
            end
        end
    end
end

