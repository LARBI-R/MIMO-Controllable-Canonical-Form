function [N] = CtrbInvMatrix(A,B)

Co = ctrb(A,B);
N=[];
n = length(A(1,:));
co = length(Co(1,:));
var = 0;

for k = 0:n-2
    for i = 0:n-2
        for j = 0:(co-n)
            Mat = Co(1+k:n-i,1+i:n+j-i-k);
            if (rank(Mat) == length(1+i:n+j-i-k) )
                N = Mat;
                var = 1;
                break;
            end
        end
        if (var ==1) 
            break;
        end
    end
    if (var ==1) 
        break;
    end
end

end
