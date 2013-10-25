function [e] = plotev(A)
    %  [e] = plotev(n)
    %
    %  computes the eigenvalues (e) of A
    %  the matrix and plots them in the complex plane.
    %

    n = size(A,1);
    e = eig(A);    % Get the eigenvalues of A

    close all    % Closes all currently open figures.
    figure(1)
    plot(real(e),imag(e),'r*') %   Plot real and imaginary parts
    xlabel('Real')
    ylabel('Imaginary')
    t1 = ['Eigenvalues of a random matrix of dimension ' num2str(n)];
    title(t1)