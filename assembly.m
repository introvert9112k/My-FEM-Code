function sctrB = assembly(e,element1)

sctr = element1(e,:);
nn   = length(sctr);

for k = 1 : nn
    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
    sctrBfem(2*k)   = 2*sctr(k)   ;
end
    sctrB = sctrBfem;
end % End of  