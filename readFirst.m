function result = readFirst(file)

% READFIRST Reads the first variable from the .mat file with the given
% name.

S = load(file);
names = fieldnames(S);
result = S.(names{1}); % Ooooh dynamic fieldnames

end

