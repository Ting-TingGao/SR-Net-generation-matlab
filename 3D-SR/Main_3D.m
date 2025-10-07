% 3D network generation main
% Author: Ting-Ting Gao
clear 
close all

alpha_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]; %thickness decay rate
D0_list = [0.05,0.2]; % initial thickness

N = 2000; % number of planes

for ai = 1:length(alpha_list)

for di = 1:length(D0_list)
for run_time = 1:5
alpha = alpha_list(ai);
D0 = D0_list(di);

num_planes = N;

vertices = [0 0 0;
            1 0 0;
            1 1 0;
            0 1 0;
            0 0 1
            1 0 1;
            1 1 1;
            0 1 1;];

k = convhull(vertices,'Simplify',true);
edges = [];
for e = 1:size(k,1)
tmp = [k(e,1),k(e,2);
k(e,2),k(e,3);
k(e,3),k(e,1);];
edges = [edges;tmp];
end
edges = removeDuplicates(edges);

faces = k;

% Initialize the cell array to store the plane coordinates
mid_plane_all = {};
upper_plane_all = {};
down_plane_all = {};
planes_coords = {};
poly_list = {};
edges_list = {};
faces_list = {};


%planes_coords{1}{1} = vertices;
poly_list{1}{1} = vertices;
edges_list{1}{1} = edges;
faces_list{1}{1} = faces;



t = 0;
trial = 1;
numP = 0;
% Finde the intersection points of the plane
while trial<1e5
    % Randomly choose part
    tmp_volume = zeros(1,length(poly_list{t+1}));
    for v = 1:length(poly_list{t+1})
        [~, volume_temp] = polyvolume(poly_list{t+1}{v});
        tmp_volume(v) = volume_temp;
    end
    prob = tmp_volume./sum(tmp_volume);
    cumulativeProbabilities = cumsum(prob);
    r = rand();
    selectedIndex = find(cumulativeProbabilities >= r, 1);
    current_vertices = poly_list{t+1}{selectedIndex};
    current_edges = edges_list{t+1}{selectedIndex};
    current_face = faces_list{t+1}{selectedIndex};
    [isValid, midplane,upper, down, planeVertices, subPoly1Vertices, subPoly2Vertices] = intersect_plane_with_polyhedron(current_vertices, current_face, D0, t, alpha);
    if isValid == 1
        mid_plane_all{t+1} = midplane;
        upper_plane_all{t+1} = upper;
        down_plane_all{t+1} = down;
        planes_coords{t+1} = planeVertices;

        
        poly_list{t+2} = poly_list{t+1};
        poly_list{t+2}(selectedIndex) = [];
        tmp_num_1 = length(poly_list{t+2});
        poly_list{t+2}{tmp_num_1+1} = subPoly1Vertices;
        poly_list{t+2}{tmp_num_1+2} = subPoly2Vertices;


        k1 = convhull(subPoly1Vertices,'Simplify',true);
        edges1 = [];
        for e = 1:size(k1,1)
        tmp = [k1(e,1),k1(e,2);
        k1(e,2),k1(e,3);
        k1(e,3),k1(e,1);];
        edges1 = [edges1;tmp];
        end
        edges1 = removeDuplicates(edges1);
        faces1 = k1;

        k2 = convhull(subPoly2Vertices,'Simplify',true);
        edges2 = [];
        for e = 1:size(k2,1)
        tmp = [k2(e,1),k2(e,2);
        k2(e,2),k2(e,3);
        k2(e,3),k2(e,1);];
        edges2 = [edges2;tmp];
        end
        edges2 = removeDuplicates(edges2);
        faces2 = k2;

        edges_list{t+2} = edges_list{t+1};
        edges_list{t+2}(selectedIndex) = [];
        tmp_num_2 = length(edges_list{t+2});
        edges_list{t+2}{tmp_num_2+1} = edges1;
        edges_list{t+2}{tmp_num_2+2} = edges2;

        faces_list{t+2} = faces_list{t+1};
        faces_list{t+2}(selectedIndex) = [];
        tmp_num_3 = length(faces_list{t+2});
        faces_list{t+2}{tmp_num_3+1} = faces1;
        faces_list{t+2}{tmp_num_3+2} = faces2;
        numP = numP + 1;
        t = t + 1;
        fprintf('Added plane number %d\n', t);
    else
        trial = trial + 1;
        %disp('Continue to search the plane');
    end

    if numP > num_planes
        disp('Done');
        break
    end
end


%% Visualization and statistical results
%visualize_planes_with_cube(planes_coords(1:100),upper_plane_all(1:100),down_plane_all(1:100));

visualize_planes_with_cube(planes_coords,upper_plane_all,down_plane_all);
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080]);
exportgraphics(gcf, 'visualized_planes_structure100.png', ...
     'BackgroundColor', 'white', 'Resolution', 300);
disp('PNG saved: visualized_planes_structure.png');
%% === Export from planes_coords ===
% fv.faces = [];
% fv.vertices = [];
% v_offset = 0;
% 
% for i = 1:length(planes_coords)
%     verts = planes_coords{i};
%     if size(verts,1) < 4
%         continue
%     end
%     try
%         k = convhull(verts, 'Simplify', true);  
%     catch
%         continue  
%     end
%     fv.vertices = [fv.vertices; verts];
%     fv.faces = [fv.faces; k + v_offset];
%     v_offset = v_offset + size(verts,1);
% end
% 
% figure;
% patch('Faces', fv.faces, 'Vertices', fv.vertices, ...
%       'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
% axis equal; view(3); camlight; lighting gouraud;
% title('Combined polyhedra from planes\_coords');

% TR = triangulation(fv.faces, fv.vertices);
% stlwrite(TR, ['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_mesh_2000.stl']);
%%
num_poly = length(planes_coords);
centroids = zeros(num_poly, 3);

for i = 1:num_poly
    centroids(i,:) = mean(planes_coords{i}, 1);
end

A = detectPolyhedronContacts(planes_coords);


%% data saving
save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_upper_planes.mat'], 'upper_plane_all');
save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_down_planes.mat'], 'down_plane_all');

save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_planes.mat'], 'planes_coords');
save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_A.mat'], 'A');
volumes_all = 0;
T = 1:1:length(planes_coords);
t = length(planes_coords);
D_matrix = D0*(T).^(-alpha);
planes_area = [];
for i = 1:t
    [~,tmpv] = polyvolume(planes_coords{i});
    area_temp = tmpv/D_matrix(i);
    planes_area(i) = area_temp;
    volumes_all = volumes_all + area_temp*D0;
end

V = 0;
for i = 1:length(planes_coords)
    [~,tmpv] = polyvolume(planes_coords{i});
V = V+tmpv;
end
save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_final_volume.mat'], 'V');
non_zero_indices = planes_area > 1e-10;
T = T(non_zero_indices);
planes_area = planes_area(non_zero_indices);
save(['.../alpha',num2str(alpha_list(ai)),'_D0_',num2str(D0_list(di)),'_',num2str(run_time),'_area.mat'], 'planes_area');


end
end

end
%% === visualizaiton ===
figure('Color','w');
num_nodes = 5;


visualize_planes_with_cube(planes_coords,upper_plane_all,down_plane_all);
hold on;

cube_vertices = [0 0 0; 1 0 0; 1 1 0; 0 1 0; ...
                 0 0 1; 1 0 1; 1 1 1; 0 1 1];
cube_edges = [1 2; 2 3; 3 4; 4 1;
              5 6; 6 7; 7 8; 8 5;
              1 5; 2 6; 3 7; 4 8];
for i = 1:size(cube_edges,1)
    v1 = cube_vertices(cube_edges(i,1), :);
    v2 = cube_vertices(cube_edges(i,2), :);
    plot3([v1(1) v2(1)], [v1(2) v2(2)], [v1(3) v2(3)], ...
        'k-', 'LineWidth', 0.02);
end
% 
% 
%[num_nodes, ~] = size(centroids);

for i = 1:num_nodes
    for j = i+1:num_nodes
        if A(i,j) == 1

            plot3([centroids(i,1), centroids(j,1)], ...
                  [centroids(i,2), centroids(j,2)], ...
                  [centroids(i,3), centroids(j,3)], ...
                  '-', 'Color', [0.9, 0.75, 0.9], 'LineWidth', 2, 'LineStyle', '-');
        end
    end
end


scatter3(centroids(1:num_nodes,1), centroids(1:num_nodes,2), centroids(1:num_nodes,3), ...
         40, 'r', 'filled');



axis equal;
axis([0 1 0 1 0 1]);
%xlabel('X'); ylabel('Y'); zlabel('Z');
axis off
view(3);
%title('Polyhedron Connectivity Network');
grid on;