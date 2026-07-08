% Author(s): Marti Geoffrey
% Epsztein Lab 2017

function wlt = fct_wavelet(x, kind, space)

switch kind
    case 'morlet'
        w0 = 6;
        switch space
            case 'time'
                t = x;
                wlt = (pi^(-0.25))*exp(1i*w0*t).*exp(-(t.^2)/2);
            case 'freq'
                wk = x;
                expnt = -(((wk - w0).^2)/2);
                wlt = (pi^(-0.25))*exp(expnt).*(wk > 0);   
        end
        
        
end




