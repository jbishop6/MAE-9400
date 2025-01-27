% % Jennifer Bishop
% % 1/26/25
% % The intention of this script is to convert .mat files to .txt files. 
% 
% clc
% clear
% close all
% 
% % Loading file
% 
% data = load('3311_surface.mat');
% 
% % Inspecting data
% 
% % disp(data)
% 
% fv = data.fv;
% % disp(data.fv)
% 
% vertices = fv.vertices;
% faces = fv.faces;
% 
% numVertices = size(vertices,1); % Number of vertices
% indices = (1:numVertices)'; % Creating indices
% 
% % Using applicable data to make same format _ref.txt file
% 
% % values = vertices(:,1);
% 
% for i = 1:numVertices
%     connected_faces = find(any(faces == i, 2));
%     values(i) = connected_faces(1);
% end
% 
% refData = [indices, values'];
% 
% 
% % size(vertices)  % Should be  [n,3]
% % size(faces) % Should be [m,3] or [m,4]
% 
% % Writing to txt file
% 
% output_file = '3311_surface.txt';
% 
% fileID = fopen(output_file,'w');
% 
% fprintf(fileID,'%6d %6d\n', refData');
% fclose(fileID)
% 
% % Assuming 'fv' is the structure containing faces and vertices
% trisurf(fv.faces, fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3));
% axis equal;




% Jennifer Bishop
% 1/26/25
% The intention of this script is to convert .mat files to .txt files. 

clc
clear
close all

% Loading file
data = load('4022_surface.mat');

% Inspecting data
fv = data.fv;
vertices = fv.vertices;
faces = fv.faces;

numVertices = size(vertices,1); % Number of vertices
indices = (1:numVertices)'; % Create indices for vertices

% Initialize values array
values = zeros(numVertices, 1);

% Loop through each vertex and assign corresponding connected face reference
for i = 1:numVertices
    % Find all faces that include vertex i
    connected_faces = find(any(faces == i, 2));  
    
    % If the vertex has at least one connected face, take the first connected face index
    if ~isempty(connected_faces)
        values(i) = connected_faces(1);  % Use the first connected face index
    else
        values(i) = 0;  % If no connection, use 0 or some placeholder
    end
end

% Combine indices and values to create the reference data
refData = [indices, values];

% Writing to .txt file
output_file = '4022_surface_ref.txt';
fileID = fopen(output_file, 'w');
fprintf(fileID, '%6d %6d\n', refData');
fclose(fileID);

% Optional: Visualize the surface mesh
trisurf(fv.faces, fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3));
axis equal;
