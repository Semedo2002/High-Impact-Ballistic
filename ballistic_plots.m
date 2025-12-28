% Analysis Plots for Ballistic Impact
% Von Mises stress, pressure, energy, and deformation plots
% Run after main simulation

clear; clc; close all;

%% Load Results
if ~isfile('ballistic_results.mat')
    error('Run ballistic_main.m first to generate results!');
end

load('ballistic_results.mat');
fprintf('Loaded results from simulation\n');

%% Extract data
time = results.time * 1000;  % Convert to milliseconds
displacement = results.displacement;
velocity = results.velocity;
stress = results.stress / 1e6;  % Convert to MPa
energy = results.energy / 1e3;  % Convert to kJ

projectile = results.projectile;
plate = results.plate;

%% Create comprehensive figure with subplots
fig = figure('Position', [50, 50, 1600, 1000], 'Color', 'w');

%% Plot 1: Displacement vs Time
subplot(3, 3, 1);
plot(time, displacement(:,3), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Z-Displacement (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('Vertical Displacement History', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 2: Velocity vs Time
subplot(3, 3, 2);
vel_mag = sqrt(sum(velocity.^2, 2));
plot(time, vel_mag, 'r-', 'LineWidth', 2);
grid on;
xlabel('Time (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Velocity (m/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impact Velocity History', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 3: Von Mises Stress vs Time
subplot(3, 3, 3);
plot(time, stress, 'Color', [0.8 0.3 0.1], 'LineWidth', 2.5);
hold on;
yline(projectile.yield/1e6, '--k', 'LineWidth', 1.5, 'Label', 'Projectile Yield');
yline(plate.yield/1e6, '--b', 'LineWidth', 1.5, 'Label', 'Plate Yield');
grid on;
xlabel('Time (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Max von Mises Stress (MPa)', 'FontSize', 11, 'FontWeight', 'bold');
title('Stress Evolution', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 10);

%% Plot 4: Kinetic Energy vs Time
subplot(3, 3, 4);
plot(time, energy, 'Color', [0.2 0.6 0.3], 'LineWidth', 2);
grid on;
xlabel('Time (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Kinetic Energy (kJ)', 'FontSize', 11, 'FontWeight', 'bold');
title('Energy Dissipation', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 5: Stress Contour (2D projection at impact time)
subplot(3, 3, 5);

% Find impact time (when velocity starts decreasing rapidly)
[~, impact_idx] = min(diff(vel_mag));
impact_time = time(impact_idx);

% Create synthetic stress field for visualization
x = linspace(-plate.length/2, plate.length/2, 50);
y = linspace(-plate.width/2, plate.width/2, 50);
[X, Y] = meshgrid(x, y);

% Radial stress distribution (higher at center)
R = sqrt(X.^2 + Y.^2);
stress_field = stress(impact_idx) * exp(-R.^2 / 0.05);

contourf(X, Y, stress_field, 20, 'LineStyle', 'none');
colormap(gca, jet);
colorbar;
axis equal;
xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Von Mises Stress at t=%.2f ms', impact_time), ...
    'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 6: Pressure Distribution
subplot(3, 3, 6);

% Calculate contact pressure (simplified)
% P = F/A, where F is impact force
contact_area = pi * projectile.radius^2;
impact_force = projectile.mass * abs(projectile.velocity) / (impact_time/1000);
pressure = impact_force / contact_area / 1e6;  % MPa

% Pressure field
pressure_field = pressure * exp(-R.^2 / 0.02);

contourf(X, Y, pressure_field, 20, 'LineStyle', 'none');
colormap(gca, hot);
colorbar;
axis equal;
xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Contact Pressure (MPa) at t=%.2f ms', impact_time), ...
    'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 7: 3D Stress Distribution
subplot(3, 3, 7);
surf(X, Y, stress_field, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap(gca, jet);
colorbar;
view(45, 30);
xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('Stress (MPa)', 'FontSize', 11, 'FontWeight', 'bold');
title('3D Stress Field', 'FontSize', 12, 'FontWeight', 'bold');
lighting gouraud;
shading interp;
set(gca, 'FontSize', 10);

%% Plot 8: Acceleration vs Time
subplot(3, 3, 8);
accel = diff(vel_mag) ./ diff(time);
accel = [accel(1); accel] / 9.81;  % Convert to g's
plot(time, accel, 'Color', [0.6 0.2 0.8], 'LineWidth', 2);
grid on;
xlabel('Time (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Deceleration (g)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impact Deceleration', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);

%% Plot 9: Damage/Penetration Analysis
subplot(3, 3, 9);

% Calculate penetration depth
penetration = -min(displacement(:,3)) * 1000;  % mm

% Create bar chart
categories = {'Max Stress', 'Plate Yield', 'Penetration'};
values = [max(stress), plate.yield/1e6, penetration*10];  % Scale penetration
colors = [0.8 0.3 0.1; 0.2 0.4 0.8; 0.3 0.7 0.3];

b = bar(values);
b.FaceColor = 'flat';
b.CData = colors;

set(gca, 'XTickLabel', categories, 'FontSize', 10);
ylabel('Value (MPa or mm×10)', 'FontSize', 11, 'FontWeight', 'bold');
title('Impact Summary', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Add value labels on bars
for i = 1:length(values)
    text(i, values(i), sprintf('  %.1f', values(i)), ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
end

%% Overall figure title
sgtitle('Ballistic Impact Analysis - High-Velocity Missile on Armored Plate', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Print Summary Statistics
fprintf('\n=== IMPACT ANALYSIS SUMMARY ===\n');
fprintf('Simulation Duration: %.3f ms\n', max(time));
fprintf('Initial Velocity: %.1f m/s\n', projectile.velocity);
fprintf('Impact Velocity: %.1f m/s\n', vel_mag(impact_idx));
fprintf('Max von Mises Stress: %.1f MPa\n', max(stress));
fprintf('Projectile Yield Strength: %.1f MPa\n', projectile.yield/1e6);
fprintf('Plate Yield Strength: %.1f MPa\n', plate.yield/1e6);
fprintf('Peak Contact Pressure: %.1f MPa\n', max(pressure_field(:)));
fprintf('Max Deceleration: %.1f g\n', max(abs(accel)));
fprintf('Estimated Penetration: %.2f mm\n', penetration);
fprintf('Initial Kinetic Energy: %.1f kJ\n', energy(1));
fprintf('Final Kinetic Energy: %.1f kJ\n', energy(end));
fprintf('Energy Absorbed: %.1f kJ (%.1f%%)\n', ...
    energy(1) - energy(end), 100*(energy(1) - energy(end))/energy(1));

% Assess damage
if max(stress) > plate.yield/1e6
    fprintf('\n*** PLATE FAILURE: Stress exceeds yield strength! ***\n');
else
    fprintf('\n*** PLATE INTACT: Stress below yield strength. ***\n');
end

fprintf('\n=== END OF ANALYSIS ===\n');