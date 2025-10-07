function [isValid, midplane, upper, down, planeVertices, subPoly1Vertices, subPoly2Vertices] = intersect_plane_with_polyhedron(vertices, faces, D0, t, alpha)
    % input:
    % vertices - polyhedron (n x 3)
    % faces - polyhedron's faces (m x k): m - number of faces，k - vertices
    % for each face
    % D0 - initial thickness
    % t - current step
    % alpha - decay exponent for thickness
    %
    % output:
    % isValid - plane is valid(0: not valid, 1: valid)
    % planeVertices - plane's vertices
    % subPoly1Vertices
    % subPoly2Vertices


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % directionType = randi([1, 2]); % 1 vertical，2 horizontal
    %  
    % if directionType == 1
    %     normal = [0, 0, 1];
    %     random_dir = randi([1, 3]);
    %     normal = circshift(normal, random_dir - 3);
    % else
    %     %normal = rand(1, 2); 
    %     normal = [1, 0, 0]; 
    %     normal = normal / norm(normal); 
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if t == 0
    %     directionType = randi([1, 2]);
    %     if directionType == 1
    %         normal = [0, 0, 1];
    %         random_dir = randi([1, 3]);
    %         normal = circshift(normal, random_dir - 3);
    %     else
    %         normal = [1, 0, 0];
    %         random_dir = randi([1, 3]);
    %         normal = circshift(normal, random_dir - 3);
    %     end
    % else
    % normal = -1 + 2 * rand(1, 3);
    % end


    normal = -1 + 2 * rand(1, 3);



    normal = normal / norm(normal); 
    min_coords = min(vertices);
    max_coords = max(vertices);

    inside = false;
    for i = 1:5000
        random_point = min_coords + (max_coords - min_coords) .* rand(1, 3);
        if inpolyhedron(faces, vertices, random_point)
            inside = true;
            point = random_point;
            break;
        end
    end

    D = D0 * (t+1)^(-alpha);
    offset1 = dot(normal, point) + D / 2;
    offset2 = dot(normal, point) - D / 2;
    mid = dot(normal, point);
    intersections1 = get_plane_intersections(vertices, faces, normal, offset1);
    intersections2 = get_plane_intersections(vertices, faces, normal, offset2);
    intersections3 = get_plane_intersections(vertices, faces, normal, mid);

    if ~isempty(intersections1) && ~isempty(intersections2) %&& size(intersections1,2) == 3 && size(intersections2,2) == 3


        [subPoly1Vertices, subPoly2Vertices] = split_polyhedron_by_planes(vertices, faces, normal, offset1, offset2);
    

        [~,originalVolume] = polyvolume(vertices);
    

        [is_v1,volume1] = polyvolume(subPoly1Vertices);
        [is_v2,volume2] = polyvolume(subPoly2Vertices);
    

        [ispoly, planeVolume] = polyvolume([intersections1;intersections2;]);
        if ispoly == 0 || is_v1 == 0 || is_v2 == 0
            isValid = 0;
            midplane = [];
            upper = [];
            down = [];
            planeVertices = [];
            subPoly1Vertices = [];
            subPoly2Vertices = []; 
        end
    
    
        planeVertices = [intersections1; intersections2];
        midplane = intersections3;
        upper = intersections1;
        down = intersections2;

    
        if abs(originalVolume - (volume1 + volume2 + planeVolume)) < 1e-7 && ispoly == 1
        %if is_v1 == 1 && is_v2 == 1 && ispoly == 1
            isValid = 1;
        else
            isValid = 0;
            midplane = [];
            upper = [];
            down = [];
            planeVertices = [];
            subPoly1Vertices = [];
            subPoly2Vertices = [];
        end
    else
        isValid = 0;
        midplane = [];
        upper = [];
        down = [];
        planeVertices = [];
        subPoly1Vertices = [];
        subPoly2Vertices = [];
    end
end

