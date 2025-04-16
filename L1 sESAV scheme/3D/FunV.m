function [ V ] = FunV( xi )

if xi > 2
    V = 0;
elseif xi <= 2 && xi >= 3/2
    V = 3/4 * ( 1 + 2 * ( 3/2 - xi )/( 3/2 - 2 ) )*( (xi - 2)/(3/2-2) ).^2 - ( xi - 3/2 ).*( (xi - 2)/(3/2-2) ).^2;
elseif xi < 3/2 && xi > 1/2
    V = xi * (2 - xi);
elseif xi <= 1/2 && xi >= 0
    V = 3/4 * ( 1 + 2 * ( (1/2-xi)/(1/2) ) ).*( ( xi - 0 )/( 1/2 ) ).^2 + ( xi - 1/2 ).*( ( xi - 0 )/( 1/2 ) ).^2;
else
    V = 0;
end


end

