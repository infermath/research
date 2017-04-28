function ur = unitrootcheck(table)
ur=0;
for i=1:size(table,1)
    if size(table,2)==1
        if table(i,1)>=1
            ur=1;
        end
    else
        % finds the roots of AR polynomial
        [x1,x2]=rootfinder(table(i,1),table(i,2));
        if (abs(x1)<=1 || abs(x2)<=1)
            ur=1; % returns 1 if the process is explosive
        end
    end
end
end