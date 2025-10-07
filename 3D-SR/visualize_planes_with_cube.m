function kk = visualize_planes_with_cube(planes, upper, down)
kk = [];
    figure;
    hold on;
    
    cube_vertices = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    cube_faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    

    num_planes = length(planes);
    blue_shades = [linspace(0.85, 0.6, num_planes)', linspace(0.85, 0.6, num_planes)', ones(num_planes, 1)];
    

    for i = 1:num_planes
        plane_vertices = planes{i};
        upper_vertices = upper{i};
        down_vertices = down{i};
        

        if size(plane_vertices, 1) < 3
            continue;
        end

        k = convhull(plane_vertices(:, 1), plane_vertices(:, 2), plane_vertices(:, 3),'Simplify',true);
        centroid_upper = mean(upper_vertices, 1);
        centroid_down = mean(down_vertices, 1);


        angles_upper = atan2(upper_vertices(:,2) - centroid_upper(2), upper_vertices(:,1) - centroid_upper(1));
        angles_down = atan2(down_vertices(:,2) - centroid_down(2), down_vertices(:,1) - centroid_down(1));
        
        [~, sort_idx_upper] = sort(angles_upper);
        [~, sort_idx_down] = sort(angles_down);
        sorted_vertices_upper = upper_vertices(sort_idx_upper, :);
        sorted_vertices_down = down_vertices(sort_idx_down, :);
        
        closed_sorted_vertices_upper = [sorted_vertices_upper; sorted_vertices_upper(1, :)];
        closed_sorted_vertices_down = [sorted_vertices_down; sorted_vertices_down(1, :)];

        faceColor = blue_shades(i, :);

        if i == 10000
            trisurf(k, plane_vertices(:, 1), plane_vertices(:, 2), plane_vertices(:, 3), ...
                'FaceColor', faceColor, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 1); 
        else
        trisurf(k, plane_vertices(:, 1), plane_vertices(:, 2), plane_vertices(:, 3), ...
                'FaceColor', faceColor, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 0); 
        end
        plot3(closed_sorted_vertices_upper(:, 1), closed_sorted_vertices_upper(:, 2), closed_sorted_vertices_upper(:, 3),'-', 'Color',[0,0,0],'LineWidth', 1);
        plot3(closed_sorted_vertices_down(:, 1), closed_sorted_vertices_down(:, 2), closed_sorted_vertices_down(:, 3), '-', 'Color',[0,0,0],'LineWidth', 1);
        kk(:,i) = length(closed_sorted_vertices_upper(:, 1));
    end
    

    patch('Vertices', cube_vertices, 'Faces', cube_faces, ...
           'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid off;
    view(3); 
    axis off

end

