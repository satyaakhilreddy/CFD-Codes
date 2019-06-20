function d=TDMA(start,iend,a,b,c,d)

    % Forward Elimination
    for i=(start+1):iend
        r=a(i)/b(i-1);
        b(i)=b(i)-r*c(i-1);
        d(i)=d(i)-r*d(i-1);
    end
    
    % Backward Substitution
    d(iend)=d(iend)/b(iend);
    for i=(iend-1):-1:start
        d(i)=(d(i)-c(i)*d(i+1))/b(i);
    end
end