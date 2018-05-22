function mesh = sampleSurfacePoints(mesh,nP)

    % input:
    % mesh - consisting of faces and vertices
    % nP    - number of points to be sampled

    mesh.area = zeros(size(mesh.faces,1),1);

    nF = size(mesh.faces,1); % number of faces
    
    %% calculate area
    for idxF = 1:nF  % for all faces

        f = mesh.faces(idxF,:);

        a = mesh.vertices( f(1),: );
        b = mesh.vertices( f(2),: );
        c = mesh.vertices( f(3),: );

        ab = b-a;
        ac = c-a;

        mesh.area(idxF) = norm(cross(ab,ac))/2;

    end
    
    %% sample surface
    % draw a face according to relative size of the faces
    p_f = mesh.area ./ sum(mesh.area);
    idxRand = discretesample( p_f, nP );
    
    % sample point on face
    r1 = rand(nP,1);
    r2 = rand(nP,1);

    A = mesh.vertices( mesh.faces(idxRand,1),: );
    B = mesh.vertices( mesh.faces(idxRand,2),: );
    C = mesh.vertices( mesh.faces(idxRand,3),: );

    mesh.points = zeros(nP,3);
    for idxP = 1:nP
      mesh.points(idxP,:) = (1 - sqrt(r1(idxP))) * A(idxP,:) + (sqrt(r1(idxP)) * (1 - r2(idxP))) * B(idxP,:) + (sqrt(r1(idxP)) * r2(idxP)) * C(idxP,:);    
    end

end
