function [acc, dt, np] =  fn_readPEER_inputFile(filename)
%% Reading from PEER format
file_in = fopen(filename, 'r');
for k = 1:3
    ans = fgetl(file_in);
end
[temp ans] = fscanf(file_in, '%c', 6);
[np ans] = fscanf(file_in, '%f', 1a);
[temp ans] = fscanf(file_in, '%c', 7);
[dt ans] = fscanf(file_in, '%f', 1);
ans = fgetl(file_in);
acc   = fscanf(file_in, '%e');
fclose(file_in);
clear k temp ans
end

