function [Ac,Bc] = ControllableFormMIMO(A,B)

P = FindMatriceP(A,B);

Ac = P * A * inv(P);
Bc = P*B;


end
