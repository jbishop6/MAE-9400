% function [pt, trg] = ReadOFF(fname, foutname)
%     % Open the file
%     fid = fopen(fname, 'r');
%     if fid == -1
%         error('Error: File %s not found.', fname);
%     end
% 
%     % Read the shape type
%     shape_type = fscanf(fid, '%s', 1);
%     if ~strcmp(shape_type, 'OFF')
%         fclose(fid);
%         error('Error: File %s is not in OFF format.', fname);
%     end
% 
%     % Read the number of points and triangles
%     a = fscanf(fid, '%f', 3);
%     if length(a) < 2
%         fclose(fid);
%         error('Error: Failed to read the number of points and triangles in %s.', fname);
%     end
%     num_pt = a(1);
%     num_trg = a(2);
% 
%     % Read point coordinates and triangle data
%     pt_coordinates = fscanf(fid, '%f', num_pt * 3);
%     temp = fscanf(fid, '%f', num_trg * 4);
% 
%     if length(pt_coordinates) ~= num_pt * 3
%         fclose(fid);
%         error('Error: Mismatch in the number of point coordinates in %s.', fname);
%     end
%     if length(temp) ~= num_trg * 4
%         fclose(fid);
%         error('Error: Mismatch in the number of triangle data in %s.', fname);
%     end
% 
%     % Reshape the data
%     pt = reshape(pt_coordinates, 3, num_pt)';
%     temp = reshape(temp, 4, num_trg)';
%     trg = temp(:, 2:4) + 1; % Adjust for 1-based indexing
% 
%     fclose(fid);
% 
%     % Optional: Write the mesh to another format and normalize it
%     if nargin == 2
%         writeObjMesh2(pt, pt, trg, foutname);
%         system(['./MeshNormal ', foutname]);
%     end
% end
function [VERT, TRIV] = ReadOFF(filename)
    % Open the file and read contents
    fid = fopen(filename, 'r');
    if fid == -1
        error('File not found');
    end
    
    % Read the header and ignore it
    header = fscanf(fid, '%s', 1);
    if ~strcmp(header, 'OFF')
        error('This is not a valid OFF file');
    end
    
    % Read the number of vertices, faces, and edges
    nverts = fscanf(fid, '%d', 1);
    nfaces = fscanf(fid, '%d', 1);
    % Skip the number of edges
    fscanf(fid, '%d', 1);

    % Read vertex coordinates
    VERT = fscanf(fid, '%f %f %f\n', [3 nverts])';

    % Read faces, expect 3 vertices per face
    TRIV = fscanf(fid, '%d %d %d\n', [3 nfaces])' + 1;  % +1 for MATLAB 1-indexing
    
    fclose(fid);
end

