function error = RMSError(V1, V2)

n = size(V1, 1);
modification = sum(V1- V2, 1) / n;
V = V2 + modification;

E = V1 - V;
E1 = E(:, 1);
E2 = E(:, 2);

error1 = sqrt(sum(E1 .* E1) / n);
error2 = sqrt(sum(E2 .* E2) / n);

error = (error1 + error2) / 2;