function dy = mechsys_01(t,y,coefficients)
%
% Input parameters:
% ------------------
% coefficients (4x1) -> a_1*d2x + a_2*d1x + a_3*d0x = 0, a_4 -> gamma
% definition (char) -> {'traditional'|'conformable'}
%
% Output parameter:
% ------------------
% ODE system ->     d1y = d2x, y(0),
%                   d1x = y,  x(0),
m = coefficients(1);
b = coefficients(2);
k = coefficients(3);

if numel(coefficients) < 4
    g = 1.0;               % If g is not included in coefficients (default)
else
    g = coefficients(4);   % Read g from the last coefficient
end
% Note: if definition = 'traditional', g is ignored

% Find auxiliar paramters like natural frequency and damping factor
omega_n = sqrt(k/m);
zeta    = b/(2*sqrt(k*m));

%From Don Juan's manuscript
        dy(1,1) = y(2);
        dy(2,1) = -2*zeta*((omega_n).^g).*(t.^(g-1)).*y(2) - ...
            ((omega_n).^(2*g)).*(t.^(2*(g-1))).*y(1);

end