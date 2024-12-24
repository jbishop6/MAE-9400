function S = read_shape_for_adj(filename)
import MESH.MESH_IO.*
fprintf('Reading mesh...'); tic;
fname = strsplit(filename,'.');
if length(fname) > 1
    if strcmp(fname{end},'obj')
        [X,T] = readObj(filename);
    elseif strcmp(fname{end},'off')
        [X,T] = readOff(filename);
    else
        error('cannot read .%s file',fname{end});
    end
else % no input file extension
    if exist([filename,'.obj'],'file')
        [X,T] = readObj(filename);
    elseif exist([filename,'.off'],'file')
        [X,T] = readOff(filename);
    else
        error('file not found: %s\n',filename)
    end
end

S.TRIV = double(T);
S.VERT = X;
S.nf = size(T,1);
S.n = size(X,1);
S.name = shape_name(filename);

S.adj = digraph(S.TRIV, S.TRIV(:, [2 3 1])); % just use matlab function seems much faster
S.adj = adjacency(S.adj);
S.adj = S.adj | S.adj';
S.adj = zeros(S.n) + S.adj;

t = toc; fprintf('done:%.4fs\n',t);
MESH.print_info(S);
end