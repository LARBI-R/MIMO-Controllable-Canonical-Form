function [N] = ComputeCob(A,B)

Ro = Indice(A,B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cob = [b_1 A^( Ro_0 - 1 )b_1 , ... , b_m A^( Ro_m - 1 )b_m ]
N=[];
m = length(B(1,:));

for i = 1:m
    
    for j = 1 : Ro(i)
        N = [N, A^(j-1) * B(:,i)];
    end
    
end

end
