
function sctrB = assembly_nonlocal(sctr,total_disp_local)

nn   = length(sctr); % Number of nodes.
sctr_n = zeros(1,nn);
for i = 1:nn
a = sctr(i);    
b = a + total_disp_local;
sctr_n(i) = b + a*2;
end

sctr_nonlocal = zeros(1,12);

cnt = 0 ;

for k = 1 : nn
    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k) - 2;

    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k) - 1;
    
    cnt = cnt + 1 ;
    sctr_nonlocal(cnt) = sctr_n(k);
end
sctrB = sctr_nonlocal;
end