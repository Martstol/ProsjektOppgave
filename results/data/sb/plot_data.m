#! /bin/octave -qf

setenv('GNUTERM', 'X11');
getenv('GNUTERM')

arg_list = argv();

if(nargin < 2)
	printf('Need name of file with data to plot and what to plot (td).\n');
	exit();
end

fid = fopen(arg_list{1}, 'r');
rawdata = fscanf(fid, '%*s %f\n');
fclose(fid);

num = length(rawdata) / 6;

printf('Number of entries: %d\n', num);

mat = vec2mat(rawdata, 6);

labels = {'advect', 'setupSolution', 'setInitialGuess', 'solve', 'project', 'windToGPU'};

advect = mat(:,1);
setupSolution = mat(:,2);
setInitialGuess = mat(:,3);
solve = mat(:,4);
project = mat(:,5);
windToGPU = mat(:,6);

data = zeros(1, 6);
data(1) = sum(advect) / num;
data(2) = sum(setupSolution) / num;
data(3) = sum(setInitialGuess) / num;
data(4) = sum(solve) / num;
data(5) = sum(project) / num;
data(6) = sum(windToGPU) / num;

exectime = sum(data);
normalizedData = data ./ exectime;

hold on;
if(strcmp(arg_list{2}, 'td'))
	pie(normalizedData, labels);
end
hold off;

print(strcat(arg_list{2}, '_', arg_list{1}, '.png'));
