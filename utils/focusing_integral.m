function [Output] = focusing_integral(f, f0, Interval, Omega)
% function [Output,Output2] = focusing_integral(f,f0,Interval,N)
% fun = @(u,n1,n2,fm) exp(1i*2*pi*u*(fm/f0*(n1-1)-(n2-1)));

for m = 1 : length(f)
    for j = 1 : size(Interval,1)
        nn1 = 1;
        for n1 = Omega'
            nn2 = 1;
            for n2 = Omega'
                w = 2*pi*(f(m)/f0*(n1-1)-(n2-1));
                if w == 0
                    Q{j}(nn1,nn2) = Interval(j,2)-Interval(j,1);
                else
                    Q{j}(nn1,nn2) = -1i/w * ( exp(1i*w*Interval(j,2)) - exp(1i*w*Interval(j,1)) ); % closed form, fast
                end
%                 Q2{j}(n1,n2) = integral(@(u)fun(u,n1,n2,f(m)),Interval(j,1),Interval(j,2));
                % for check
            nn2 = nn2 + 1;
            end
            nn1 = nn1 + 1;
        end
    end
    Output{m} = sum(cat(3,Q{:}),3);
%         Output2{m} = sum(cat(3,Q2{:}),3);
end

end