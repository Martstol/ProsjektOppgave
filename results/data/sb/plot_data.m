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
plot(cuda(1, :), cuda(2, :), "Color", "blue");
plot(petsc_cpu_cg(1, :), petsc_cpu_cg(2, :), "Color", "green");
plot(petsc_cpu_gmres(1, :), petsc_cpu_gmres(2, :), "Color", "red");
plot(petsc_gpu_cg(1, :), petsc_gpu_cg(2, :), "Color", "magenta");
plot(petsc_gpu_gmres(1, :), petsc_gpu_gmres(2, :), "Color", "black");
xlabel("Grid points", "FontSize", 12, "FontName", "Arial");
set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
ylabel("Execution time (seconds)", "FontSize", 12, "FontName", "Arial");
legend("Cuda (old)", "PETSc CPU CG", "PETSc CPU GMRES", "PETSc GPU CG", "PETSc GPU GMRES",
	"Location", "SouthOutside");
hold off;
print("exec_time_all.png");

clf;

hold on;
title("GPU Solver Execution Time", "FontSize", 12, "FontName", "Arial");
plot(cuda(1, :), cuda(2, :), "Color", "blue");
plot(petsc_gpu_cg(1, :), petsc_gpu_cg(2, :), "Color", "magenta");
plot(petsc_gpu_gmres(1, :), petsc_gpu_gmres(2, :), "Color", "black");
xlabel("Grid points", "FontSize", 12, "FontName", "Arial");
set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
ylabel("Execution time (seconds)", "FontSize", 12, "FontName", "Arial");
legend("Cuda (old)", "PETSc GPU CG", "PETSc GPU GMRES", "Location", "SouthOutside");
hold off;
print("exec_time_gpu.png");

