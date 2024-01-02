    % Define Earth properties
    earth_radius = 6371;  % Earth's radius in kilometers
    
    % Create a circle representing Earth
    earth_theta = linspace(0, 2*pi, 1000);  % Angle range for plotting
    earth_x = earth_radius * cos(earth_theta);
    earth_y = earth_radius * sin(earth_theta);
    
    % Plot Earth centered at (0,0)
    figure;
    plot(earth_x, earth_y, 'r');  % Plot Earth as a red circle centered at (0,0)
    hold on;  % Hold the plot for adding the GTO, launch point, and rocket trajectory
    
    % Define the orbital elements for GTO
    semi_major_axis = 24322.0;  % Example value for semi-major axis in kilometers
    eccentricity = 0.7;          % Example value for eccentricity
    
    % Calculate the parameters needed for plotting the GTO ellipse
    theta = linspace(0, 2*pi, 1000);  % Angle range for plotting
    r = semi_major_axis * (1 - eccentricity^2) ./ (1 + eccentricity * cos(theta)); % Polar equation for ellipse
    
    % Convert polar coordinates to Cartesian coordinates for GTO
    x_gto = r .* cos(theta);
    y_gto = r .* sin(theta);
    
    % Plot the GTO ellipse centered at (-1.7025e+04,0)
    plot(x_gto, y_gto,'g');
    xlabel('X-axis (km)');
    ylabel('Y-axis (km)');
    title('Earth, Geostationary Transfer Orbit, and Rocket Trajectory from Kourou');
    legend('Earth', 'GTO Orbit');
    
    % Add the marker for Kourou launch site
    kourou_latitude = 5.2360;   % Latitude of Kourou
    kourou_longitude = -52.7686;  % Longitude of Kourou
    
    % Convert Kourou coordinates to Cartesian coordinates on the Earth's surface
    x_kourou = earth_radius * cosd(kourou_latitude) * cosd(kourou_longitude);
    y_kourou = earth_radius * cosd(kourou_latitude) * sind(kourou_longitude);
    
    % Plotting Kourou on the Earth's surface
    plot(x_kourou, y_kourou, 'ro', 'MarkerSize', 10);  % Plotting Kourou as a red circle marker
    text(x_kourou, y_kourou, 'Kourou', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'); % Label for Kourou
    
    % Calculate the angle along the GTO orbit for the rocket's destination
    %rocket_destination_angle = pi/4; % Example angle, in radians (modify as needed)
    
    % Calculate the coordinates of the rocket's destination on the GTO orbit
    %x_rocket_dest = semi_major_axis * cos(rocket_destination_angle);
    %y_rocket_dest = semi_major_axis * sin(rocket_destination_angle);
    
    % Plotting the trajectory from Kourou to the rocket's destination
    %plot([x_kourou, x_rocket_dest], [y_kourou, y_rocket_dest], 'b--', 'LineWidth', 1.5); % Green dashed line for rocket trajectory
    
    axis equal; % Set aspect ratio to equal for a circular representation of the ellipse
    % (Your previous code remains unchanged)
    
    % Define rocket parameters
    rocket_height = 5000; % Example height reached by the rocket (in kilometers)
    rocket_distance = 3500; % Example distance traveled horizontally (in kilometers)
    
    % Calculate the rocket's initial and final positions
    x_rocket_start = x_kourou;
    y_rocket_start = y_kourou;
    

    % Define the parabolic trajectory coordinates starting from the rocket's initial position
    x_parabola = linspace(0, rocket_distance, 100);
    y_parabola = (rocket_height / (rocket_distance^2)) * x_parabola.^2;

    % Offset the parabolic trajectory to match the rocket's starting point
    x_parabola = x_parabola + x_rocket_start;
    y_parabola = y_parabola + y_rocket_start;

  
    % Plotting the rocket's parabolic trajectory
    plot(x_parabola, y_parabola, 'b--', 'LineWidth', 1.5); % Blue dashed line for rocket's parabolic trajectory


    grid on;
    hold off;  % Release the plot hold
