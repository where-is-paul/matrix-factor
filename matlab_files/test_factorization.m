%tests the factorization generated by matrix-factor
clc;
mat_names = [];
for i = 1:8
    mat_names = [mat_names; strcat('testmat', num2str(i))];
end
mat_names = cellstr(mat_names);
other_mats = { 'aug3dcqp_rcm_perm'; 'bloweya'; 'bratu3d'; 'tuma1'; ...
            '1138_bus'};

all_mats = other_mats;%mat_names;%[mat_names; other_mats];
        
lfil = 2;
tol = 0.001;
for i = 1:length(all_mats)
    mat_name = all_mats{i};
    fprintf('Now testing %s:\n', mat_name);
    base = '';
    file = strcat(base, mat_name, '.mtx');
    %base = 'C:\Users\Paul\Documents\My Dropbox\Online Resources\My
    %Homework\UBC 2011W\NSERC\matrix-factor-recode\';

    cmd = horzcat('ldl_driver ', num2str(lfil), ' ', num2str(tol));
    cmd = horzcat(cmd, ' test_matrices/', file, ' -y -n');
    [~, ~] = system(cmd);

    %pause(0.5);

    A = mmread(file);
    l = mmread(strcat(base, 'outL.mtx'));
    d = mmread(strcat(base, 'outD.mtx'));
    p = mmread(strcat(base, 'outPerm.mtx'));
    p = diag(p)+1;

    %avg_col_norm = norm(A,1)/size(A,1);
    %the residual only provides a heuristic. e.g. lfil = 5, droptol = 0.1 for
    %testmat 5 gives 4.5*10^-2 residual but lfil = 5, droptol = 0.01 gives
    %1.8*10^-1 as a residual. the latter is a better preconditioner, but
    %measures further from the original matrix.
    fprintf('The relative residual is %d.\n', norm(A(p,p) - l*d*l', 1)/norm(A, 1));
    %spy(A(p,p)); figure; spy(abs(l*d*l') > 1e-8);

    e = ones(size(A,1),1);
    gmres(A(p,p),e,60,1e-6,250,l*d, l');
    %gmres(A,e,30,1e-6,2);
    fprintf('\n');
end
