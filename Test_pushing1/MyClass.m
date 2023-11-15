classdef MyClass
   % Some attributes are set to logical values
   properties (Constant = true)
      CnstProp = 2*pi
      Prop2 = MyClass.setupAccount
   end
   properties
      % Static method of this class
      Prop1 = MyClass.setupAccount
      x = MyClass.test1(2)
      %[rx,ry]= MyClass.robot_shape
      % Constant property from this class
%       Prop2 = MyClass.CnstProp
%       % Function that returns a value
%       Prop3 = datestr(now)
%       % A class constructor
%       Prop4 = AccountManager
   end
   methods (Static)
      function accNum = setupAccount
         accNum = randi(9,[1,12]);
      end
      function [r_x, r_y]=robot_shape
                     % khối xe
                     r = 0.1; 
                     t_s = 0:pi/180:2*pi;
                     r_x = r*cos(t_s);
                     r_y  = r*sin(t_s);
      end

   end
   methods (Static)
      function x = test1(a)
                x = 0;
                [y1,y2] = test2(3)
               x = y1 +y2;
          function [x,y] = test2(b)
              y = 2+b;
              x = 1+b;
          end
      end
      

   end
   
   
end