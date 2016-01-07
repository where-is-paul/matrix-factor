function x = run_ildl( L, D, p, S, b )
  %
  % x = run_ildl( L, D, p, S, b )
  %
  % return an approximate solution x of the system Ax = b
  % given the incomplete factorization P'SASP = LDL'.
  %

  Sb = S * b;
  PSb = Sb(p);
  x1 = L' \ (D \ (L \ PSb));
  x1(p) = x1;
  x = S * x1;
end

