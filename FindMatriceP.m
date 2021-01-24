function [P] = FindMatriceP(A,B)

P =[];


m = length(B(1,:));
n = length(A(1,:));

Sigma = zeros(1,n);

Cob = ComputeCob(A,B);
Ro = Indice(A,B);

for i = 1:m
    for j = 1:i
        Sigma(i) = Sigma(i) + Ro(j);
    end
end

CobInv = inv(Cob);

for i = 1:m
    for j = 1: Ro(i)
        P = [P;CobInv(Sigma(i),:)*A^(j-1)];
    end
end

end
