% 3D Dynamic Visualization of Ballistic Impact
% Animates projectile impact with color-coded stress
% Run after main simulation

clear; clc; close all;

%% Load Results
if ~isfile('ballistic_results.mat')
    error('Run ballistic_main.m first to generate results!');
end

load('ballistic_results.mat');
fprintf('Loaded results from simulation\n');

%% Reconstruct mesh for visualization
projectile = results.projectile;
plate = results.plate;
nodes = results.nodes;
n_proj_nodes = results.n_proj_nodes;
n_total_nodes = size(nodes, 1);

% Time steps to animate
n_frames = min(100, length(results.time));
frame_indices = round(linspace(1, length(results.time), n_frames));

%% Setup Figure
fig = figure('Position', [100, 100, 1400, 800], 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

% Set view angle
view(ax, 45, 30);
axis(ax, 'equal');
xlabel(ax, 'X (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'Y (m)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel(ax, 'Z (m)', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, 'Ballistic Missile Impact on Armored Plate', 'FontSize', 14, 'FontWeight', 'bold');

% Lighting
light('Position', [1 1 1], 'Style', 'infinite');
light('Position', [-1 -1 1], 'Style', 'infinite');
lighting(ax, 'gouraud');
material(ax, 'metal');

% Colormap for stress
colormap(ax, jet);
c = colorbar(ax);
c.Label.String = 'Von Mises Stress (MPa)';
c.FontSize = 11;

%% geometry

ntheta = 20;
nz = 20;
theta = linspace(0, 2*pi, ntheta);
z_cyl = linspace(0, projectile.length, nz);
[Theta, Z] = meshgrid(theta, z_cyl);
X_proj = projectile.radius * cos(Theta);
Y_proj = projectile.radius * sin(Theta);
Z_proj = Z;

% plate viz
x_plate = linspace(-plate.length/2, plate.length/2, 30);
y_plate = linspace(-plate.width/2, plate.width/2, 30);
[X_plate, Y_plate] = meshgrid(x_plate, y_plate);
Z_plate_top = (projectile.length + 0.1) * ones(size(X_plate));
Z_plate_bot = Z_plate_top - plate.thickness;

%surface objects
h_proj = surf(ax, X_proj, Y_proj, Z_proj);
h_plate_top = surf(ax, X_plate, Y_plate, Z_plate_top);
h_plate_bot = surf(ax, X_plate, Y_plate, Z_plate_bot);

set(h_proj, 'EdgeColor', 'none', 'FaceColor', 'interp');
set(h_plate_top, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.9);
set(h_plate_bot, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.7);

x_ground = linspace(-plate.length, plate.length, 10);
y_ground = linspace(-plate.width, plate.width, 10);
[X_ground, Y_ground] = meshgrid(x_ground, y_ground);
Z_ground = (projectile.length + 0.1 + plate.thickness + 0.05) * ones(size(X_ground));
h_ground = surf(ax, X_ground, Y_ground, Z_ground);
set(h_ground, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Text
txt_time = text(ax, -0.4, -0.4, projectile.length + 0.5, '', ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'Margin', 5);

txt_vel = text(ax, -0.4, -0.4, projectile.length + 0.45, '', ...
    'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

txt_stress = text(ax, -0.4, -0.4, projectile.length + 0.4, '', ...
    'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

%% Animation
fprintf('Starting animation...\n');

for frame = 1:n_frames
    idx = frame_indices(frame);
    t = results.time(idx);
    progress = idx / length(results.time);
    z_disp_proj = projectile.velocity * t;
    Z_proj_new = Z_proj + z_disp_proj;
    
    % stress
    if results.stress(idx) > 0
        stress_val = min(results.stress(idx) / 1e6, 1000); % MPa, capped at 1000
    else
        stress_val = 0;
    end
    
    % stress color
    impact_factor = max(0, 1 - abs(Z_proj_new(:) - (projectile.length + 0.1)) / 0.2);
    stress_colors = stress_val * impact_factor;
    stress_colors = reshape(stress_colors, size(Z_proj_new));
    
    % Update projectile
    set(h_proj, 'XData', X_proj, 'YData', Y_proj, 'ZData', Z_proj_new, ...
        'CData', stress_colors);
    
    % exagg. deformation 
    if progress > 0.3  
        contact_factor = (progress - 0.3) / 0.7;
        
        % crater
        R_plate = sqrt(X_plate.^2 + Y_plate.^2);
        deform = -0.02 * contact_factor * exp(-R_plate.^2 / 0.05);
        
        Z_plate_top_def = Z_plate_top + deform;
        Z_plate_bot_def = Z_plate_bot + 0.5 * deform;
        
        % plate stress
        plate_stress = stress_val * exp(-R_plate.^2 / 0.1);
        
        set(h_plate_top, 'ZData', Z_plate_top_def, 'CData', plate_stress);
        set(h_plate_bot, 'ZData', Z_plate_bot_def);
    end
    
    % Update text
    set(txt_time, 'String', sprintf('Time: %.3f ms', t*1000));
    set(txt_vel, 'String', sprintf('Velocity: %.1f m/s', ...
        projectile.velocity + 9.81*t));
    set(txt_stress, 'String', sprintf('Max Stress: %.0f MPa', stress_val));
    
    caxis(ax, [0, max(500, stress_val)]);
    
    if progress < 0.5
        zlim(ax, [-0.1, projectile.length + 0.6]);
    else
        zlim(ax, [projectile.length - 0.1, projectile.length + 0.4]);
    end
    
    xlim(ax, [-0.6, 0.6]);
    ylim(ax, [-0.6, 0.6]);
    
    drawnow;
    
    % Write frame to video (if enabled)
    % frame_img = getframe(fig);
    % writeVideo(v, frame_img);
    
    pause(0.05);  % Adjust for animation speed
end


fprintf('Animation complete!\n');
fprintf('Run ballistic_plots.m for detailed analysis plots\n');