function G = Godunov( ul,ur )
%The Godunov function
U=1-(ul+ur);
%Case 1
if ul<0.5&&ur<0.5
    us=ul;
%Case 2
else if ul>0.5&&ur>0.5
        us=ur;
%Case 4
else if ul>0.5&&ur<0.5
        us=1/2;
%Case 3
else if U<0
        us=ul;
else
        us=ur;
            end
        end
    end
end

G=us*(1-us);
