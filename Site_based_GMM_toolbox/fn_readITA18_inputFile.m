function [acc, dt, np] =  fn_readITA18_inputFile(filename)
% Open the file
file_in = fopen(filename, 'r');


% dt
while true
    line = fgetl(file_in);
    if contains(line, 'SAMPLING_INTERVAL_S:')
        break;  % Found the USER5 line
    end
end
% Extract sampling interval after ':'
parts = split(line, ':');
dt = str2double(parts{2});

% acc
while true
    line = fgetl(file_in);
    if contains(line, 'USER5:')
        break;  % Found the USER5 line
    end
end
ans = fgetl(file_in);
acc   = fscanf(file_in, '%e');
np = numel(acc);

fclose(file_in);
end