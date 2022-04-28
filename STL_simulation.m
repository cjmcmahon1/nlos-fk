data = stlread("3DBenchy.stl");
obj_coords = downsample(data.Points, 10);
xvals = linspace(-120, 0, 32);
yvals = ones(1, 32) .* 60;
zvals = linspace(-40, 80, 32);
[X, Z] = meshgrid(xvals, zvals);
figure(1)
trimesh(data, 'FaceColor', 'none', 'EdgeColor', 'k');
xlabel('X'); ylabel('Y'); zlabel('Z');
hold on;
wall_coords = [X(:) 60.*ones(size(X(:), 1), 1) Z(:)];
% scatter3(X(:), 60.*ones(size(X(:), 1), 1), Z(:), 'b');
scatter3(wall_coords(:,1), wall_coords(:,2), wall_coords(:,3), 'b');
diffx = obj_coords(1, 1) - xvals(16);
diffy = obj_coords(1, 2) - yvals(16);
diffz = obj_coords(1, 3) - zvals(16);
lines = line(wall_coords, pi/4);
linedist = distance_to_line(obj_coords(1,:), lines(:,:,1), lines(:,:,2));

threshold = 1;
plotall_lines = false;
for l = 1:size(lines, 1)
    obj_distances = distance_to_line(obj_coords, ...
                                     lines(l, :, 1), lines(l, :, 2));
    obj_coords_on_line = obj_coords(obj_distances < threshold, :);
    if size(obj_coords_on_line, 1) >= 1
        %plot the vector of the shortest distance ray that hits the object
        vec_to_wall = obj_coords_on_line - lines(l, :, 1);
        [min_val, min_idx] = min(vecnorm(vec_to_wall, 2, 2));
        end_coord = obj_coords_on_line(min_idx,:);
        pltline(lines(l, :, 1), end_coord, 0.8)
    elseif plotall_lines
        pltline(lines(l, :, 1), lines(l, :, 2), 0.1)
    end
end

hold off;

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