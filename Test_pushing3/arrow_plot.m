function poly = arrow_plot(pos , head, scale)
vertices = [-0.5, 0; -0.5, 2; -1.5, 2; 0, 3; 1.5,2; 0.5,2; 0.5,0]'*scale;


matrix_t = [cos(head- pi/2) -sin(head -pi/2) pos(1) ;...
              sin(head-pi/2) cos(head -pi/2) pos(2);...
              0            0        1   ];
          ver = [];
 for i = 1:length(vertices)
    no = (matrix_t)*[vertices(1,i); vertices(2,i);1];
    ver(1,i) = no(1);
    ver(2,i) = no(2);
 end

% Create a polyshape
poly = polyshape(ver');

% Plot the polyshape

plot(poly);

end
