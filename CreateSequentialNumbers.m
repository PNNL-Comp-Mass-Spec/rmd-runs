function [Int_Vector] = CreateSequentialNumbers(Int_Vector)





a = unique(Int_Vector);
while(length(a) < max(Int_Vector))            
    k = 1;
    while(a(k) == k)
        k = k+1;
    end
    Int_Vector(Int_Vector == a(k)) = k;
    a = unique(Int_Vector);
end
