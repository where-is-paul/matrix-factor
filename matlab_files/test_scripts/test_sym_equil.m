%tests the factorization generated by matrix-factor
disp('======================== test start ======================');
mat_names = [];
for i = 1:8
    mat_names = [mat_names; strcat('testmat', num2str(i))];
end
mat_names = cellstr(mat_names);
other_mats = { 'aug3dcqp'; 'bloweya'; 'bratu3d'; ...
            'stokes64'; 'tuma1'; 'tuma2'; ...
            '1138_bus'};

all_mats = mat_names;%{'tuma2'};%[mat_names; other_mats];%
        
lfil = 3;
tol = 0.00;
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

    A = mmread(file);
    B = mmread(strcat(base, 'outB.mtx')); 
    S = mmread(strcat(base, 'outS.mtx'));
    p = mmread(strcat(base, 'outPerm.mtx'));
    p = diag(p);
    
    P = speye(size(A));
    P = P(:,p);
   
    fprintf('The residual between scaled and original is: %d\n', ...
        norm(B - P'*S*A*S*P, 1)/norm(B,1));

    bool_equil = true;
    for j = 1:size(B,1)
        if (abs(norm(B(:,j), Inf)-1) > 1e-8) && (abs(norm(B(j,:), Inf)-1) > 1e-8)
            bool_equil = false;
            break;
        end
    end

    if (~bool_equil)
        fprintf('Matrix not equilibriated.\n');
        fprintf('Check index %d.\n', j);
    else
        fprintf('Matrix equilibriated.\n');
    end
    fprintf('The condition number is %d.\n', condest(B));
    fprintf('\n');
end

disp('======================== test end ======================');