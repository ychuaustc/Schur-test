function error = RMSError1(U1, U2)

n = size(U1, 2);

modification = sum(U1- U2) / n;
U = U2 + modification;
E = U1 - U;

% E = U1 - U2;

error = sqrt(sum(E .* E) / n);