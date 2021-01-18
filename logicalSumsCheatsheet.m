%% Summing up logical arrays cheatsheet

sum(logicals) > 0 % True if ANY are TRUE

sum(logicals) == 0 % True if ALL are FALSE

sum(logicals) < length(logicals) % True if ANY are FALSE

sum(logicals) == length(logicals) % True if ALL are TRUE