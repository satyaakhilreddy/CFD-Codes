%% Function - minmod

function s=minmod(a,b)
    if a*b>0
        s=(a/abs(a))*min(abs(a),abs(b));
    else
        s=0;
    end
end