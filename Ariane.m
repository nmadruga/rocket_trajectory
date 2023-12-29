%function rocket_simulation()
    % System parameters
    m = 1000;          % Rocket mass (kg)
    g = 9.81;          % Gravity acceleration (m/s^2)
    mu = 3.986e14;     % Standard gravitational parameter of Earth (m^3/s^2)
    I = 5000;          % Rocket moment of inertia (kg*m^2)
    % Initial conditions
    V0 = 0;            % Initial velocity (m/s)
    a_gamma = deg2rad(90);       % Initial launch angle (degrees to radiant)
    y0 = 0;            % Initial altitude (m)
    x0 = 0;            % Initial downrange distance (m)
    M0 = m;            % Initial mass (kg)
    atvc =deg2rad(2);% TVC angle
    % Integration options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    % Integrate the differential equations
    [t, state] = ode45(@rocket_eqns, [0, 100], [V0, a_gamma, y0, x0, M0], options, m, g, mu, I, atvc)
        % Extract results
    V = state(:, 1);
    a_gamma = (state(:, 2));
    y = state(:, 3);
    x = state(:, 4);
    M = state(:, 5);
    % Display results
    figure;
    subplot(2, 2, 1);
    plot(t, V);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Rocket Velocity');
    subplot(2, 2, 2);
    plot(t, a_gamma);
    xlabel('Time (s)');
    ylabel('Launch Angle (degrees)');
    title('Rocket Launch Angle');
    subplot(2, 2, 3);
    plot(t, y);
    xlabel('Time (s)');
    ylabel('Altitude (m)');
    title('Rocket Altitude');
    subplot(2, 2, 4);
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Downrange Distance (m)');
    title('Rocket Downrange Distance');
    % Differential equations function
    function dydt = rocket_eqns(t, state, m, g, mu, I, atvc);
        V = state(1);
        a_gamma = state(2);
        y = state(3);
        x = state(4);
        M = state(5);
        r = sqrt(x^2 + y^2);
        T = thrust_force(t);   % Function for thrust force as a function of time
        D = drag_force(V);     % Function for drag force as a function of velocity
        L = lift_force(V);     % Function for lift force as a function of velocity
        dydt = zeros(5, 1);
        dydt(1) = (T*cos(atvc) - D)/m  - mu*sin(a_gamma)/(r^2) 
        dydt(2) = (V*cos(a_gamma)/r + T*sin(atvc)/(m*V) + L/(m*V) - mu*cos(a_gamma)/(V*r^2))
        dydt(3) = V*sin(a_gamma)
        dydt(4) = V*cos(a_gamma)
        dydt(5) = -T/(I*g)
    end
    % Function for thrust force as a function of time
    function T = thrust_force(t)
        % Example function based on time
        T = 500*sin(2*pi*0.1*t);
    end
    % Function for drag force as a function of velocity
    function D = drag_force(V)
        % Specific parameters for your design
        rho = 1.2;       % Air density (kg/m^3)
        CD = 0.2;        % Drag coefficient
        A = 5;           % Reference cross-sectional area of the rocket (m^2)
        % Drag force formula
        D = 0.5 * rho * V^2 * CD * A;
    end
    % Function for lift force as a function of velocity
    function L = lift_force(V)
        % Specific parameters for your design
        rho = 1.2;       % Air density (kg/m^3)
        CL = 0.1;        % Lift coefficient
        A = 5;           % Reference cross-sectional area of the rocket (m^2)
        % Lift force formula
        L = 0.5 * rho * V^2 * CL * A;
    end
%end