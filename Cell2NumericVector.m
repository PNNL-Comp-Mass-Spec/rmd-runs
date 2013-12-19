function [Int_Vector] = Cell2NumericVector(Cell_Vector)

n = length(Cell_Vector);
[b,c] = unique(Cell_Vector);
m = length(b);
a = b;
c = sort(c);
for i = 1:m
    a(i) = Cell_Vector(c(i));
end
Int_Vector = zeros(n,1);
for i = 1:m
    b = strmatch(a(i),Cell_Vector);
    Int_Vector(b) = i;
end

