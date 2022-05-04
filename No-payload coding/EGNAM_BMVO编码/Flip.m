function [simplexOut,fval,func_evals] = Flip(funfcn,simplex,fv,dim,lb,ub)
LB = ones(dim,1)*lb;
UB = ones(dim,1)*ub;

n = 30;
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;



%v = zeros(n,n+1);
%fv = zeros(1,n+1);


%v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
%x = xin;    % Change x to the form expected by funfcn
%fv(:,1) = funfcn(x,varargin{:});
%fv(:,1) = funfcn(x);

func_evals = 0;
itercount = 0;
how = '';


v = simplex;

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
%usual_delta = 0.05;             % 5 percent deltas for non-zero terms
%zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x


%for j = 1:n
%    y = xin;
%    if y(j) ~= 0
%        y(j) = (1 + usual_delta)*y(j);
%    else
%        y(j) = zero_term_delta;
%    end
%    v(:,j+1) = y;
%    x(:) = y; f = funfcn(x);
%    fv(1,j+1) = f;
%end

% sort so v(1,:) has the lowest function value

[fv,j] = sort(fv);

v = simplex(:, j);

 xbar = sum(v(:,one2n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    %ADDED BOUND
    xr = ApplyBound(xr,LB, UB);

    x = xr;

    %fxr = funfcn(x,varargin{:});
    fxr = funfcn(x);
    
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xbar=xbar(:);
        et = (1 + rho*chi)*xbar;
        ek = rho*chi*v(:,end);
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
    
       %ADDED BOUND
        xe = ApplyBound(xe, LB, UB);

        x(:) = xe;

    
        %        fxe = funfcn(x,varargin{:});
        fxe = funfcn(x);
        
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);

                %ADDED BOUND
                xc = ApplyBound(xc, LB, UB);

                x(:) = xc;
                
%                fxc = funfcn(x,varargin{:});
                fxc = funfcn(x);

                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                %ADDED BOUND
                xcc = ApplyBound(xcc, LB, UB);

                
                x(:) = xcc; fxcc = funfcn(x);
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=two2np1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    
                    x(:) = v(:,j); 
                     %ADDED BOUND
                   x = ApplyBound(x, LB, UB);
                    v(:,j) = x; 
                    
                    fv(:,j) = funfcn(x);
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fval,j] = sort(fv);
    simplexOut = v(:,j);
   

end 