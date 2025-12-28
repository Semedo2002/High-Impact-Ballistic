% FEA Main Ballistic Impact FEA Simulation
clear; clc; close all;
%% Params
% Projectile
projectile.length = 0.3;
projectile.radius = 0.025;
projectile.mass = 5.0;
projectile.velocity = -800;
projectile.E = 200e9;
projectile.nu = 0.3;% Poisson
projectile.rho = 7850;
projectile.yield = 350e6;

%plate
plate.length = 1.0;
plate.width = 1.0;
plate.thickness = 0.05;
plate.E = 210e9;% Steel
plate.nu = 0.3;
plate.rho = 7850;
plate.yield = 500e6;
plate.mass = plate.rho * plate.length * plate.width * plate.thickness;  % Calculate mass

%% Mesh
fprintf('Generating mesh...\n');

% Projectile
nz_proj = 30;
nr_proj = 8;
ntheta = 16; %tip ref

% Plate
nx_plate = 40;
ny_plate = 40;
nz_plate = 10;%ref at impct

%nodes
proj_nodes = [];
proj_elements = [];
node_id = 1;

%projectile
for i = 1:nz_proj
    z = (i-1) * projectile.length / (nz_proj-1);
    for j = 1:nr_proj
        r = (j-1) * projectile.radius / (nr_proj-1);
        for k = 1:ntheta
            theta = (k-1) * 2*pi / ntheta;
            x = r * cos(theta);
            y = r * sin(theta);
            proj_nodes(node_id, :) = [x, y, z];
            node_id = node_id + 1;
        end
    end
end

%nodes
plate_nodes = [];
z_offset = projectile.length + 0.1; % 10cm gap initially

for i = 1:nx_plate
    x = (i-1) * plate.length / (nx_plate-1) - plate.length/2;
    for j = 1:ny_plate
        y = (j-1) * plate.width / (ny_plate-1) - plate.width/2;
        for k = 1:nz_plate
            z = z_offset + (k-1) * plate.thickness / (nz_plate-1);
            plate_nodes(node_id, :) = [x, y, z];
            node_id = node_id + 1;
        end
    end
end
all_nodes = [proj_nodes; plate_nodes];
n_proj_nodes = size(proj_nodes, 1);
n_plate_nodes = size(plate_nodes, 1);
n_total_nodes = n_proj_nodes + n_plate_nodes;

fprintf('Total nodes: %d\n', n_total_nodes);
fprintf('Projectile nodes: %d, Plate nodes: %d\n', n_proj_nodes, n_plate_nodes);

%% State Vectors
u = zeros(n_total_nodes, 3);
v = zeros(n_total_nodes, 3);
a = zeros(n_total_nodes, 3);
v(1:n_proj_nodes, 3) = projectile.velocity;

%% Params
dt = 1e-6;
t_total = 0.003;                 % 3 millisecond simulation
n_steps = round(t_total / dt);
save_interval = 50;

% Storage for results
time_history = [];
displacement_history = [];
velocity_history = [];
stress_history = [];
energy_history = [];

%% TI
fprintf('Starting simulation...\n');
fprintf('Total steps: %d\n', n_steps);
% Newmark
beta = 0.25;
gamma = 0.5;
% Contact params
contact_stiffness = 1e10;
contact_damping = 1e6;

for step = 1:n_steps
    t = step * dt;
    
    % forces
    F_ext = zeros(n_total_nodes, 3);
    F_gravity = zeros(n_total_nodes, 3);
    F_gravity(1:n_proj_nodes, 3) = -projectile.mass * 9.81 / n_proj_nodes;
    F_gravity(n_proj_nodes+1:end, 3) = -plate.mass * 9.81 / n_plate_nodes;
    
    % forces penalty 
    F_contact = zeros(n_total_nodes, 3);
    
    for i = 1:n_proj_nodes
        proj_pos = all_nodes(i, :) + u(i, :);
        
        % contact check?
        for j = n_proj_nodes+1:n_total_nodes
            plate_pos = all_nodes(j, :) + u(j, :);
            
            dist_vec = proj_pos - plate_pos;
            dist = norm(dist_vec);
            
            if dist < 0.01  
                
                normal = dist_vec / (dist + eps);
                
                penetration = max(0, 0.01 - dist);
                
                % Penalty force
                F_n = contact_stiffness * penetration;
                
                % Damping force
                rel_vel = v(i, :) - v(j, :);
                F_d = contact_damping * dot(rel_vel, normal);
                
                F_total = (F_n + F_d) * normal;
                
                F_contact(i, :) = F_contact(i, :) - F_total;
                F_contact(j, :) = F_contact(j, :) + F_total;
            end
        end
    end

    F_ext = F_gravity + F_contact;
    
    % lumped mass 
    M = zeros(n_total_nodes, 1);
    M(1:n_proj_nodes) = projectile.mass / n_proj_nodes;
    M(n_proj_nodes+1:end) = plate.rho * (plate.length * plate.width * plate.thickness) / n_plate_nodes;

    for i = 1:n_total_nodes
        a(i, :) = F_ext(i, :) / M(i);
    end

    plate_bottom_idx = n_proj_nodes + (1:nx_plate*ny_plate);
    a(plate_bottom_idx, :) = 0;
    v(plate_bottom_idx, :) = 0;
    
    % Newmark integr
    u = u + dt * v + (0.5 - beta) * dt^2 * a;
    v_old = v;
    v = v + (1 - gamma) * dt * a;
    
    all_nodes_current = all_nodes + u;
    
    v = v + gamma * dt * a;
    
    % Save results at intervals
    if mod(step, save_interval) == 0
        time_history(end+1) = t;
        displacement_history(end+1, :) = [mean(u(:,1)), mean(u(:,2)), mean(u(:,3))];
        velocity_history(end+1, :) = [mean(v(:,1)), mean(v(:,2)), mean(v(:,3))];
        
        % von misses stress
        stress = calculate_stress(u, all_nodes, projectile, plate, n_proj_nodes);
        stress_history(end+1) = max(stress);
        
        KE = 0.5 * sum(M .* sum(v.^2, 2));
        energy_history(end+1) = KE;
        
        if mod(step, save_interval*10) == 0
            fprintf('Step %d/%d (%.1f%%), t=%.4f ms, Max stress=%.1f MPa\n', ...
                step, n_steps, 100*step/n_steps, t*1000, max(stress)/1e6);
        end
    end
end

%% Results
results.time = time_history;
results.displacement = displacement_history;
results.velocity = velocity_history;
results.stress = stress_history;
results.energy = energy_history;
results.nodes = all_nodes;
results.u_final = u;
results.projectile = projectile;
results.plate = plate;
results.n_proj_nodes = n_proj_nodes;

save('ballistic_results.mat', 'results');
fprintf('Results saved to ballistic_results.mat\n');

fprintf('\n=== Simulation Complete ===\n');
fprintf('Run visualization scripts:\n');
fprintf('  - ballistic_visualize.m for 3D animation\n');
fprintf('  - ballistic_plots.m for analysis plots\n');

%% Stress
function stress = calculate_stress(u, nodes, projectile, plate, n_proj_nodes)
    
    n_nodes = size(nodes, 1);
    stress = zeros(n_nodes, 1);
    
    for i = 1:n_nodes

        if i <= n_proj_nodes
            E = projectile.E;
        else
            E = plate.E;
        end
        % Characteristic length
        strain = norm(u(i, :)) / 0.01;
        stress(i) = E * strain;
    end
end