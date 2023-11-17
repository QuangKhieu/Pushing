classdef Inner_control < handle
    %INNER_CONTROL Summary of this class goes here
    % inner control to maining pushing, exert desinged force f_id
    %   Detailed explanation goes here
    
    properties
        K_omega = 0.10;
        K_Pin = 0.04;
        K_Iin = 0.05;
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
                    continue;       
                end
                
                e_fi = f_di(i) - f_i ;
                
                %PI control
                robs(i).omega = obj.K_omega*robs(i).theta;
                v     = obj.K_Pin*e_fi + ...
                                obj.K_Iin*(obj.old_error(i) + e_fi*dt);
                %chuan hoa  
                if v > 0.05
                     robs(i).v = 0.05;
                elseif(v<0)
                    robs(i).v = 0.01;
                else
                    robs(i).v = v;
                end
                    
            end
        end
    end
end

