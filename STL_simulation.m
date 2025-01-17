%%
data = stlread("3DBenchy.stl");
obj_coords = downsample(data.Points, 10);
%reorder object coordinates by z-height for reflectance
% obj_coords = sort(obj_coords, 3); %not actually needed
%organize reflectance such that objects at increasing z height have
%increasing reflectance
reflectances = zeros(size(obj_coords, 1), 1);
reflectances(obj_coords(:,3) < 10) = 1;
reflectances((obj_coords(:,3) < 30) & (obj_coords(:,3) >= 10)) = 0.2;
reflectances((obj_coords(:,3) >= 30)) = 1;
obj_coords = cat(2, obj_coords, reflectances);
%generate coordinates of wall
wall_res = 64;
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
%generate parallel lines coming from each wall coordinate
lines = line(wall_coords, pi/4);

threshold = 3; %max distance between ray & object point to consider "hit"
timing_res = 256; %number of bins in timing histogram
max_ray_length = 5; %parameter to scale the ray length into a time 
plotall_lines = false;
time_coords = zeros(wall_res^2, timing_res);
ray_distances = zeros(wall_res^2, 1);
c = 1.5e8; %speed of light halved to account for the doubled path
for l = 1:size(lines, 1)
    %for each line, find all object points sufficiently close to that line
    %(within the above threshold). Then, pick the point with the shortest
    %euclidean distance to the wall as the reflection point. 
    obj_distances = distance_to_line(obj_coords(:,1:3), ...
                                     lines(l, :, 1), lines(l, :, 2));
    obj_coords_on_line = obj_coords(obj_distances < threshold, :);
    %make a histogram of timing values, 1 if a hit was detected, 0 if no
    %hits were detected
    timing_hist = zeros(timing_res, 1); 
    if size(obj_coords_on_line, 1) >= 1
        %plot the vector of the shortest distance ray that hits the object
        vec_to_wall = obj_coords_on_line(:, 1:3) - lines(l, :, 1);
        [min_val, min_idx] = min(vecnorm(vec_to_wall, 2, 2));
        %store the (x, y, t) measurement using the ray distance and the
        %speed of light
        end_coord = obj_coords_on_line(min_idx,:);
        reflectance = end_coord(4);
        %set the hit to be at the closest histogram bin
        rel_distance = (sqrt(min_val)-7) / max_ray_length;
        timing_hist(round(rel_distance * timing_res)) = reflectance;
        pltline(lines(l, :, 1), end_coord(1:3), 0.8)
        ray_distances(l) = sqrt(min_val);
    elseif plotall_lines
        pltline(lines(l, :, 1), lines(l, :, 2), 0.1)
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

function line_coords = line(start_coords, xy_angle)
    l = 150;
    s = start_coords;
    e = start_coords + [l*cos(xy_angle) -1*l*sin(xy_angle) 0];
    line_coords = cat(3, s, e);
end

function d = distance_to_line(pt, line_start, line_end)
    v1 = pt - line_start;
    v2 = pt - line_end;
    denom = vecnorm(line_end - line_start, 2, 2);
    num = vecnorm(cross(v1, v2, 2), 2, 2);
    d = num/denom;
end