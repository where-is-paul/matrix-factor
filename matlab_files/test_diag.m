i=1;
while (i < size(A,1))
    if (abs(d(i,i)) < 1e-12 && abs(d(i, i+1)) < 1e-12) 
        if(abs(d(i, i+1)) < 1e-12 && min(abs(d(i,i)), abs(d(i+1,i+1))) < 1e-12)
            disp(i);
            break;
        end
    end
    
    if (abs(d(i,i+1)) ~= 0)
        i = i+2;
    else
        i = i+1;
    end
end

for i = 1:size(A,1)
    if norm(L(:,i) - l(:,i), 1) > 1e-6
        disp(i);
        disp(norm(L(:,i) - l(:,i), 1));
        break;
    end
end