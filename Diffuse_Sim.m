% Load an STL file and extract the vertices as a point cloud
data = stlread("3DBenchy.stl");
obj_coords = downsample(data.Points, 5);
%reorder object coordinates by z-height for reflectance
% obj_coords = sort(obj_coords, 3); %not actually needed
%organize reflectance such that objects at increasing z height have
%increasing reflectance
reflectances = zeros(size(obj_coords, 1), 1);
reflectances(obj_coords(:,3) < 10) = 10;
reflectances((obj_coords(:,3) < 30) & (obj_coords(:,3) >= 10)) = 10;
reflectances((obj_coords(:,3) >= 30)) = 10;
obj_coords = cat(2, obj_coords, reflectances);
%generate coordinates of wall
wall_res = 32;
xvals = linspace(-120, 0, wall_res);
yvals = ones(1, 32) .* 60;
zvals = linspace(-40, 80, wall_res);
[X, Z] = meshgrid(xvals, zvals);
figure(1)
trimesh(data, 'FaceColor', 'none', 'EdgeColor', 'k');
xlabel('X'); ylabel('Y'); zlabel('Z');
hold on;
wall_coords = [X(:) 60.*ones(size(X(:), 1), 1) Z(:)];
scatter3(wall_coords(:,1), wall_coords(:,2), wall_coords(:,3), 'b');

threshold = 1; %max distance between ray & object point to consider "hit"
timing_res = 256; %number of bins in timing histogram
max_ray_length = 15; %parameter to scale the ray length into a time 
offset = 7;
plotall_lines = false;
time_coords = zeros(wall_res^2, timing_res);
c = 1.5e8; %speed of light halved to account for the doubled path
for l = 1:size(wall_coords, 1)
    %make a histogram of timing values, 1 if a hit was detected, 0 if no
    %hits were detected
    timing_hist = zeros(timing_res, 1); 
    %for each line, find all object points sufficiently close to that line
    %(within the above threshold). Then, pick the point with the shortest
    %euclidean distance to the wall as the reflection point.
    %generate rays deflecting at varying angles off of the wall
    rays = diffuserays(wall_coords(l,:), pi/4);
    ray_distances = zeros(size(rays, 1), 1);
    for r = 1:size(rays, 1)
        %for each ray, find distance to closest point, and then convert
        %that distance into a timing measurement, scaled by the reflectance
        %of the point hit. We scan over multiple rays at different angles
        %and increment the histogram bins for each ray, for each point on
        %the wall.
        obj_distances = distance_to_line(obj_coords(:,1:3), ...
                                         rays(r, :, 1), rays(r, :, 2));
        obj_coords_on_line = obj_coords(obj_distances < threshold, :);
        if size(obj_coords_on_line, 1) >= 1
            %plot the vector of the shortest distance ray that hits the object
            vec_to_wall = obj_coords_on_line(:, 1:3) - rays(r, :, 1);
            [min_val, min_idx] = min(vecnorm(vec_to_wall, 2, 2));
            %store the (x, y, t) measurement using the ray distance and the
            %speed of light
            end_coord = obj_coords_on_line(min_idx,:);
            reflectance = end_coord(4);
            %set the hit to be at the closest histogram bin
            rel_distance = (sqrt(min_val)-offset) / max_ray_length;
            hist_idx = round(rel_distance * timing_res);
            timing_hist(hist_idx) = timing_hist(hist_idx) + reflectance;
            %pltline(rays(r, :, 1), end_coord(1:3), 0.1)
            ray_distances(r) = sqrt(min_val);
        end
    end
    time_coords(l, :) = timing_hist;
end
meas_sim = flip(reshape(time_coords, wall_res, wall_res, timing_res), 1);
hold off;

% run f-k migration
wall_size = 2; % scanned area is 2 m x 2 m
tof = ones(wall_res, wall_res); %tof is arbitrarily set to one in perfect sim
fprintf('\nRunning f-k migration\n');
fk = cnlos_reconstruction(meas_sim, tof, wall_size, 2, 256);
fbp = cnlos_reconstruction(meas_sim, tof, wall_size, 0, 256);
lct = cnlos_reconstruction(meas_sim, tof, wall_size, 1, 256);
function pltline(s, e, alpha)
    plot3([s(1), e(1)], [s(2) e(2)], [s(3) e(3)], ...
        'LineWidth', 0.1, 'Color', [1 0 0 alpha])
end

function line_coords = diffuserays(start_coord, xy_angle)
    l = 150;
    d_theta = pi/40;
    num_angles = 10;
    tvals = linspace(xy_angle-d_theta, xy_angle+d_theta, num_angles);
    pvals = linspace(-d_theta, d_theta, num_angles);
    [THETA, PHI] = meshgrid(tvals, pvals);
    s = ones(size(THETA(:), 1), 3) .* start_coord;
    addition = [l*cos(THETA(:)) -1*l*sin(THETA(:)) l*sin(PHI(:))];
    e = start_coord + addition;
    line_coords = cat(3, s, e);
end

function d = distance_to_line(pt, line_start, line_end)
    v1 = pt - line_start;
    v2 = pt - line_end;
    denom = vecnorm(line_end - line_start, 2, 2);
    num = vecnorm(cross(v1, v2, 2), 2, 2);
    d = num/denom;
end