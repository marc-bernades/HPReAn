function    F = Filter_Matrix(N)


%% Second order filter (Reference Lele 1992 JCP)
alpha = 0.4627507;
beta  = 0.2265509;

a = 1/16*(11 + 10*alpha - 10*beta);
b = 1/32*(15 + 34*alpha + 30*beta);
c = 1/16*(-3 + 6*alpha  + 26*beta);
d = 1/32*(1  - 2*alpha  + 2*beta);

% Inatialize matricies each direction
A = zeros(N,N);
B = zeros(N,N);

% Set filter matrices for each direction

% Sweep internal points
for j = 4:N-3
    A(j,j) = a;
    A(j,j+3) = d/2;
    A(j,j-3) = d/2;
    A(j,j+2) = c/2;
    A(j,j-2) = c/2;
    A(j,j+1) = b/2;
    A(j,j-1) = b/2;

    B(j,j)   = 1;
    B(j,j+2) = beta;
    B(j,j-2) = beta;
    B(j,j+1) = alpha;
    B(j,j-1) = alpha;
end

% Boundaries
A(1,1) = 15/16;
A(1,2) = 1/16*4;
A(1,3) = 1/16*-6;
A(1,4) = 1/16*4;
A(1,5) = 1/16*-1;

A(2,1) = 1/16;
A(2,2) = 3/4;
A(2,3) = 1/16*6;
A(2,4) = 1/16*-4;
A(2,5) = 1/16*1;

A(3,1) = -1/16;
A(3,2) = 1/16*4;
A(3,3) = 5/8;
A(3,4) = 1/16*4;
A(3,5) = 1/16*-1;

A(end,1) = 15/16;
A(end,2) = 1/16*4;
A(end,3) = 1/16*-6;
A(end,4) = 1/16*4;
A(end,5) = 1/16*-1;

A(end-1,1) = 1/16;
A(end-1,2) = 3/4;
A(end-1,3) = 1/16*6;
A(end-1,4) = 1/16*-4;
A(end-1,5) = 1/16*1;

A(end-2,1) = -1/16;
A(end-2,2) = 1/16*4;
A(end-2,3) = 5/8;
A(end-2,4) = 1/16*4;
A(end-2,5) = 1/16*-1;

B(1,1) = 1; B(2,2) = 1; B(3,3) = 1;
B(end,end) = 1; B(end-1,end-1) = 1; B(end-2,end-2) = 1;

F = inv(B)*A;


end