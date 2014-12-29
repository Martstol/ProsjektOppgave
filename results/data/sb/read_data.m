function [points, data] = read_data(filename)

fid = fopen(filename, 'r');
rawdata = fscanf(fid, '%*s %f\n');
fclose(fid);

points = rawdata(1);
rawdata(1) = [];

data = sum(rawdata) / length(rawdata);
