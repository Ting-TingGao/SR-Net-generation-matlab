function [ispolyhedron, V] = polyvolume(vertices)
    vertices = unique_with_threshold(vertices, 1e-10);

    if size(vertices, 1) < 4
        ispolyhedron = 0;
        V = NaN;  
        return;  
    else
        ispolyhedron = 1;
    end

    if rank(vertices(2:end, :) - vertices(1, :)) < 3
        ispolyhedron = 0;
        V = NaN;
        disp('The points are coplanar or collinear.');
        return;
    end
    

    K = convhull(vertices(:,1), vertices(:,2), vertices(:,3),'Simplify',true);
    

    V = 0;
    

    ref_point = mean(vertices);
    

    for i = 1:size(K,1)

        face_vertices = vertices(K(i,:), :);
        

        tetrahedron = [ref_point; face_vertices];
        

        a = tetrahedron(2,:) - tetrahedron(1,:);
        b = tetrahedron(3,:) - tetrahedron(1,:);
        c = tetrahedron(4,:) - tetrahedron(1,:);
        tetra_volume = abs(dot(a, cross(b,c))) / 6;
        
        V = V + tetra_volume;
    end
end

