function sqrtA = square_root(A, eps)
    C = chol(A);
    [u, h, its] = poldec(C, eps);
    sqrtA = h;
end