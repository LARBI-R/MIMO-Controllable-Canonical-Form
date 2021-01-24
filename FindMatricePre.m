function [R] = FindMatricePre(A,B)

R =[];

P = FindMatriceP(A,B);
Bc = P * B;
Ro = Indice(A,B);
m = length(B(1,:));

Bcc = zeros(sum(Ro),m);

somme = 0;
somme2 = Ro(1);

Bcc(1:Ro(1),1) = Bc(1:Ro(1),1);

for i = 2:m
    somme2 = somme2+Ro(i);
    somme = somme+Ro(i-1);
    Bcc(1+somme : somme2,i) =  Bc(1+somme :  somme2,i);
end

R = inv( (Bc' * Bc) ) * Bc' * Bcc;

end
