% Pure roll error estimation
a_y = linspace(0.1, 3, 100) * 9.81;             % lateral acceleration (m/s^2)
roll_gradient = (linspace(0.5, 2, 100) / 9.81)';   % deg roll per m/s^2

roll_angle = roll_gradient * a_y;                % in degrees

% What the IMU reads
a_y_imu = a_y .* cos(deg2rad(roll_angle));

% Error of the IMU in percentage
a_y_error = (a_y - a_y_imu) ./ a_y * 100;

figure;
surf(roll_gradient, a_y, a_y_error)
title('Ay Error vs Roll gradient vs Lateral Acceleration')
xlabel('Roll Gradient Deg per m/s^2')
ylabel('Lateral Acceleration m/s^2')
zlabel('IMU acceleration Error %')
