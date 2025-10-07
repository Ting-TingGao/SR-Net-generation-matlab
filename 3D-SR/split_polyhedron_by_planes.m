function [subPoly1Vertices, subPoly2Vertices] = split_polyhedron_by_planes(vertices, faces, normal, offset1, offset2)
    intersections1 = get_plane_intersections(vertices, faces, normal, offset1);
    intersections2 = get_plane_intersections(vertices, faces, normal, offset2);

    subPoly1Vertices = [];
    subPoly2Vertices = [];

    for i = 1:size(vertices, 1)
        if dot(normal, vertices(i, :) - (normal * offset1)) > 0
            subPoly1Vertices = [subPoly1Vertices; vertices(i, :)];
        else
            subPoly2Vertices = [subPoly2Vertices; vertices(i, :)];
        end
    end
    
    subPoly1Vertices = [subPoly1Vertices; intersections1];
    subPoly2Vertices = [subPoly2Vertices; intersections2];
end