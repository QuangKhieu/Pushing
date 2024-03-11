% Parameters for the ellipse
center = [0.5, 0.5];  % Center of the ellipse
axesLengths = [0.4, 0.2];  % Major and minor axes lengths
rotationAngle = 30;  % Rotation angle in degrees

% Calculate ellipse vertices
theta = linspace(0, 2*pi, 100);
x = center(1) + axesLengths(1) * cos(theta);
y = center(2) + axesLengths(2) * sin(theta);

% Rotate ellipse vertices
R = [cosd(rotationAngle), -sind(rotationAngle); sind(rotationAngle), cosd(rotationAngle)];
ellipseVertices = R * [x; y];

% Plot the rotated ellipse and fill with color using the fill function
figure;
obs = Obs_elipse();
tic
fill (obs.so_xy(1,:), obs.so_xy(2,:), 'red');
toc
figure(2)
tic
fill(ellipseVertices(1, :), ellipseVertices(2, :), 'cyan');
toc
axis equal;
title('Filled Ellipse Shape using rectangle function');
