function [h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=selectCoefficients(method,params,varargin)
    methodName=method;
    switch method
        case 'exact1'
            if length(params)~=2
                error('Wrong number of parameters: needs 2, %d given',...
                    length(params))
            end
            s = 1;
            phi = @(x,v) exp(-(x.^2+v.^2)./(2.*s));%./(2.*pi.*s);
            h = @(x,v) 0;
            fx = @(x,v) -v;
            fv = @(x,v) 0;
            gxx = @(x,v) 0;
            gvv = @(x,v) params(1);
            gxv = @(x,v) 0;
            sx = @(x,v) 0;
            sv = @(x,v) params(2);
            hasExact = true;
        case 'Langevin1'
            s = 1;
            phi = @(x,v) exp(-(x.^2+v.^2)./(2.*s));%./(2.*pi.*s);
            h = @(x,v) 0;
            fx = @(x,v) -v;
            fv = @(x,v) 0;
            gxx = @(x,v) 0;
            gvv = @(x,v) params(1).*(1+1./(1+x.^2));
            gxv = @(x,v) 0;
            sx = @(x,v) 0;
            sv = @(x,v) params(2).*sqrt(1+1./(1+x.^2));
            hasExact = false;
        otherwise
            error('Unknown coefficients: %s given',method)
    end
end