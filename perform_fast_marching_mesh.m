function [D, S, Q] = perform_fast_marching_mesh(vertex, faces, start_points, options)
    % perform_fast_marching_mesh - launch the Fast Marching algorithm on a 3D mesh.
    %
    %   [D, S, Q] = perform_fast_marching_mesh(vertex, faces, start_points, options)
    %
    %   vertex, faces: a 3D mesh
    %   start_points(i) is the index of the ith starting point.
    %
    %   D is the distance function to the set of starting points.
    %   S is the final state of the points: -1 for dead, 0 for open, 1 for far.
    %   Q is the index of the closest point. Q is set to 0 for far points.
    %       Q provides a Voronoi decomposition of the domain.
    %
    %   Optional:
    %   - You can provide non-uniform speed in options.W.
    %   - You can provide special conditions for stop in options:
    %       'options.end_points': stop when these points are reached.
    %       'options.nb_iter_max': stop when a given number of iterations is reached.
    %   - You can provide a heuristic in options.heuristic (typically that tries to guess the distance
    %       that remains from a given node to a given target).
    %       This is an array of the same size as W.
    %   - You can provide a map L=options.constraint_map that reduces the set of
    %       explored points. Only points with current distance smaller than L
    %       will be expanded. Set some entries of L to -Inf to avoid any
    %       exploration of these points.
    %
    %   Copyright (c) 2004-2006 Gabriel Peyrè

    % Initialize options
    options.null = 0;
    nverts = max(size(vertex));
    
    end_points = getoptions(options, 'end_points', []);
    verbose = getoptions(options, 'verbose', 1);
    nb_iter_max = getoptions(options, 'nb_iter_max', Inf);
    W = getoptions(options, 'W', ones(nverts, 1));
    L = getoptions(options, 'constraint_map', []);
    H = getoptions(options, 'heuristic', []);
    values = getoptions(options, 'values', []);
    dmax = getoptions(options, 'dmax', 1e8);
    
    I = find(L == -Inf); L(I) = -1e9;
    I = find(L == Inf); L(I) = 1e9;
    
    nb_iter_max = min(nb_iter_max, 1.2 * max(size(W)));
    
    if size(vertex, 1) > size(vertex, 2)
        vertex = vertex';
    end
    if size(faces, 1) > size(faces, 2)
        faces = faces';
    end
    start_points = start_points(:);
    end_points = end_points(:);
    
    % Display input sizes for debugging
    fprintf('vertex size: %d x %d\n', size(vertex, 1), size(vertex, 2));
    fprintf('faces size: %d x %d\n', size(faces, 1), size(faces, 2));
    fprintf('start_points size: %d x %d\n', size(start_points, 1), size(start_points, 2));
    fprintf('end_points size: %d x %d\n', size(end_points, 1), size(end_points, 2));
    fprintf('W size: %d x %d\n', size(W, 1), size(W, 2));
    
    % Call the fast C-coded version if possible
    try
        [D, S, Q] = perform_front_propagation_mesh(vertex, faces - 1, W, start_points - 1, end_points - 1, nb_iter_max, H, L, values, dmax);
        Q = Q + 1;
    catch ME
        disp('Error occurred in perform_front_propagation_mesh:');
        disp(ME.message);
        rethrow(ME);
    end
    
    % Replace C 'Inf' value (1e9) by Matlab Inf value.
    D(D > 1e8) = 0;
    D = sparse(D);
    D = reshape(D, length(vertex), length(start_points));
end
