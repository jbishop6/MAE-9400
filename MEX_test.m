% Define example input data with correctly transposed vertices and faces
vertex = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]; % Transpose to 3 x nverts
faces = [1, 1, 1; 2, 3, 2; 3, 4, 4]; % Correct transposition to 3 x nfaces
start_points = [1];
nverts = size(vertex, 2); % Number of vertices

% Ensure vertex is correctly transposed to 3 x nverts
if size(vertex, 1) > size(vertex, 2)
    vertex = vertex';
end

% Ensure faces is correctly transposed to 3 x nfaces
if size(faces, 1) > size(faces, 2)
    faces = faces';
end

% Adjust options to match transposed vertex
options.W = ones(nverts, 1); % Non-uniform speed array
options.end_points = [];
options.nb_iter_max = 30;
options.constraint_map = ones(nverts, 1); % Constraint map
options.heuristic = [];
options.values = [];
options.dmax = 1e8;

% Display inputs for verification
disp('vertex:'); disp(vertex);
disp('faces:'); disp(faces);
disp('start_points:'); disp(start_points);
disp('options:'); disp(options);

% Call the MEX file
try
    [D, S, Q] = perform_front_propagation_mesh(vertex, faces - 1, options.W, start_points - 1, options.end_points - 1, options.nb_iter_max, options.heuristic, options.constraint_map, options.values, options.dmax);
    disp('Distance matrix D:'); disp(D);
    disp('State matrix S:'); disp(S);
    disp('Closest points Q:'); disp(Q);
catch ME
    disp('An error occurred:');
    disp(ME.message);
end
