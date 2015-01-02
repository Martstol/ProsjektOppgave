#! /bin/octave -qf

setenv('GNUTERM', 'X11');
getenv('GNUTERM')

n = 7;

cuda = zeros(2, n);
petsc_cpu = zeros(2, n);
petsc_gpu = zeros(2, n);

for i = 1:n
	[p, d] = read_data(strcat("cuda_conf", int2str(i-1)));
	cuda(1, i) = p;
	cuda(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_cpu_conf", int2str(i-1)));
	petsc_cpu(1, i) = p;
	petsc_cpu(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_gpu_conf", int2str(i-1)));
	petsc_gpu(1, i) = p;
	petsc_gpu(2, i) = d;
end

hold on;
title("Wind Simulation Execution Time", "FontSize", 12, "FontName", "Arial");

loglog(cuda(1, :), cuda(2, :), "Color", "blue", "LineWidth", 1);
loglog(petsc_cpu(1, :), petsc_cpu(2, :), "Color", "green", "LineWidth", 1);
loglog(petsc_gpu(1, :), petsc_gpu(2, :), "Color", "red", "LineWidth", 1);

xlabel("Grid points", "FontSize", 12, "FontName", "Arial");
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

ylabel("Execution time (seconds)", "FontSize", 12, "FontName", "Arial");
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));

legend("SOR (CUDA)", "PETSc CPU", "PETSc GPU", "Location", "SouthOutside");

hold off;
print("exec_time_all.png");

