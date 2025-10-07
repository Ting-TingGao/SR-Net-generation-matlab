function intersections = get_plane_intersections(vertices, faces, normal, offset)
    intersections = [];
    
    % Define a helper function to check if a point is inside a triangle
    function inside = is_point_in_triangle(pt, v1, v2, v3)
        % Barycentric coordinates method
        u = v2 - v1;
        v = v3 - v1;
        w = pt - v1;
        uv = dot(u, v);
        uu = dot(u, u);
        vv = dot(v, v);
        wu = dot(w, u);
        wv = dot(w, v);
        D = uv * uv - uu * vv;
        if D == 0
            inside = false;
            return;
        end
        s = (uv * wv - vv * wu) / D;
        tt = (uv * wu - uu * wv) / D;
        inside = (s >= 0) && (tt >= 0) && (s + tt <= 1);
    end

    % Loop over each face (triangle) and check intersections
    for i = 1:size(faces, 1)
        face = faces(i, :);
        v1 = vertices(face(1), :);
        v2 = vertices(face(2), :);
        v3 = vertices(face(3), :);

        % Calculate the intersection points with each edge of the triangle
        edges = [v1; v2; v2; v3; v3; v1];
        for j = 1:size(edges, 1) / 2
            v1_edge = edges(2 * j - 1, :);
            v2_edge = edges(2 * j, :);

            % Compute the intersection of edge with the plane
            t = (offset - dot(normal, v1_edge)) / dot(normal, v2_edge - v1_edge);
            if t >= 0 && t <= 1
                intersection = v1_edge + t * (v2_edge - v1_edge);

                % Check if intersection is inside the triangle
                if is_point_in_triangle(intersection, v1, v2, v3)
                    intersections = [intersections; intersection];
                end
            end
        end
    end
end