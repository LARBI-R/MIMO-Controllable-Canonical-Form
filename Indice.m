function [ c ] = Indice(A,B)

C = ctrb(A,B);

n = length(A(1,:));
m = length(B(1,:));
co = length(C(1,:));

c = zeros(1,m);

ind = [];

mtest = 0;
cas2 = 0;
degTest = 0;
degTest2 = 0;

for i = 0 : co-m
        Mattest = C(1:n,1:1+m+i);

        if ( rank(Mattest) == n )
            for k = 1: n-1
                for f = 1:m
                    if ( C(:,1+m+i) == A^(k)*B(:,f) )
                    mtest = f;
                    degTest = k;
                    break;
                    end
                end
            if (mtest ~= 0)    
                break;
            end
            end
        end
        
        if ( rank(Mattest) < n )
            for k = 1: n-1
                for f = 1:m
                    if ( C(:,1+m+i) == A^(k)*B(:,f) )
                    degTest2 = k;
                    ind = [ind, f];
                    cas2 = 1;
                    break;
                    end
                end
             if (cas2 == 1)    
                break;
             end
            end
        end
        
        if(mtest ~= 0) 
            break;
        end
end

   

c(mtest) = degTest+1;


if ( mtest > 1)

    for i = 1:mtest-1
        c(i) = degTest+1;
    end

    for i = mtest+1:m
        c(i) = degTest;
    end

end

if ( mtest == 1)

    for i = mtest+1:m
        c(i) = degTest;
    end
    
end

if ( cas2 == 1)
    for i = 1:length(ind)
        c(ind(i)) = degTest2;
    end
end

end