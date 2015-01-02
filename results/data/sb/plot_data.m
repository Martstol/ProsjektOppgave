#! /bin/octave -qf

setenv('GNUTERM', 'X11');
getenv('GNUTERM')

n = 7;

cuda = zeros(2, n);
petsc_cpu_cg = zeros(2, n);
petsc_cpu_gmres = zeros(2, n);
petsc_gpu_cg = zeros(2, n);
petsc_gpu_gmres = zeros(2, n);

for i = 1:n
	[p, d] = read_data(strcat("cuda_", int2str(i-1)));
	cuda(1, i) = p;
	cuda(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_cpu_cg_", int2str(i-1)));
	petsc_cpu_cg(1, i) = p;
	petsc_cpu_cg(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_cpu_gmres_", int2str(i-1)));
	petsc_cpu_gmres(1, i) = p;
	petsc_cpu_gmres(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_gpu_cg_", int2str(i-1)));
	petsc_gpu_cg(1, i) = p;
	petsc_gpu_cg(2, i) = d;
	
	[p, d] = read_data(strcat("petsc_gpu_gmres_", int2str(i-1)));
	petsc_gpu_gmres(1, i) = p;
	petsc_gpu_gmres(2, i) = d;
end

hold on;
title("Solver Execution Time", "FontSize", 12, "FontName", "Arial");

loglog(cuda(1, :), cuda(2, :), "Color", "black", "LineWidth", 1, "-.");
loglog(petsc_cpu_cg(1, :), petsc_cpu_cg(2, :), "Color", "green", "LineWidth", 1, "-+");
loglog(petsc_cpu_gmres(1, :), petsc_cpu_gmres(2, :), "Color", "magenta", "LineWidth", 1, "-*");
loglog(petsc_gpu_cg(1, :), petsc_gpu_cg(2, :), "Color", "red", "LineWidth", 1, "-o");
loglog(petsc_gpu_gmres(1, :), petsc_gpu_gmres(2, :), "Color", "blue", "LineWidth", 1, "-x");

xlabel("Grid points", "FontSize", 12, "FontName", "Arial");
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

ylabel("Execution time (seconds)", "FontSize", 12, "FontName", "Arial");
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));

legend("SOR (CUDA)", "PETSc CPU CG", "PETSc CPU GMRES", "PETSc GPU CG", "PETSc GPU GMRES", "Location", "SouthOutside");

hold off;
print("exec_time_all.png");

clf;

speedup_cg = zeros(2, n);
speedup_gmres = zeros(2, n);

for i = 1:n
	speedup_cg(1, i) = petsc_gpu_cg(1, i);
	speedup_cg(2, i) = petsc_cpu_cg(2, i) / petsc_gpu_cg(2, i);
	
	speedup_gmres(1, i) = petsc_gpu_gmres(1, i);
	speedup_gmres(2, i) = petsc_cpu_gmres(2, i) / petsc_gpu_gmres(2, i);
end

hold on;
title("Solver Speedup", "FontSize", 12, "FontName", "Arial");

semilogx(speedup_cg(1, :), speedup_cg(2, :), "Color", "blue", "LineWidth", 1);
semilogx(speedup_gmres(1, :), speedup_gmres(2, :), "Color", "red", "LineWidth", 1);

xlabel("Grid points", "FontSize", 12, "FontName", "Arial");
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

ylabel("Speedup", "FontSize", 12, "FontName", "Arial");
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));

legend("CG", "GMRES", "Location", "SouthOutside");

hold off;
print("speedup_all.png");

