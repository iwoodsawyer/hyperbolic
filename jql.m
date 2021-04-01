function [Q,R,J] = jql(A,J,econ,tol)
%JQL J-orthogonal (or J-unitary, or hyperbolic) QL factorization
%   [Q,L] = JQL(A,J) computes a J-orthogonal (or J-unitary, or hyperbolic)
%   QL factorization of A such that matrix A = Q*L, where L is an lower
%   triangular matrix with positive values on the diagonal, and Q is a
%   J-orthogonal matrix, where Q'*J*Q = J. The signature matrix J must be a
%   matrix with 1 or -1 on the diagonal and zeros on the subdiagonals.
%
%   [Q,L,Jp] = JLQ(A,J) returns a permuted J-orthogonal, which is the case
%   when row permutations are applied for the hyperbolic rotations. If row
%   permutations are applied, than Q'*J*Q = Jp.
%
%   [Q,L,Jp] = JQL(A,J,0) produces the "economy size" decomposition.
%
%   [Q,L,Jp] = JQL(A,J,[],tol) is an optional tolerance parameter can be
%   specified. It is used as an absolute threshold to detect zero elements
%   during the reduction process. The default tol is REALMAX/EPS.
%
%   The algorithm consists in zeroing iteratively entries in each column of
%   A by application of orthogonal (plane) rotations, row permutations and
%   J-orthogonal (hyperbolic) rotations, see [1]. The algorithm is updated
%   to include the forward stable hyperbolic rotations in [2].
%
%   References:
%   [1] D. Henrion, P. Hippe, Hyperbolic QR factorization for J-spectral
%   factorization of polynomial matrices. LAAS-CNRS Research Report,
%   Toulouse, France, February 2003.
%
%   [2] S. Chandrasekaran, A.H. Sayed, Stabilizing the generalized Schur
%   algorithm, SIAM J. Matrix Anal. Appl., Vol. 17, No. 4, pp. 950-983,
%   October 1996.
%
%   [3] E. Anderson, Discontinuous Plane Rotations and the Symmetric
%   Eigenvalue Problem, Technical Report, December, 2000.
%
%   See also: QR.

% Based on: http://homepages.laas.fr/henrion/software/hqr/hqr.m
% which is written by D. Henrion and P. Hippe, November 22, 2002.

% Revised by: Ivo Houtzager [2011]
% -	better numerical implementation of hyperbolic and Givens rotations, see [2] [3].
% - better selection of tolerance
% - support for economy size


% check the number of input arguments
if nargin < 1
    error('JQL requires at least two input arguments')
end

% assign default values to unspecified parameters
[m,n] = size(A);
if (nargin < 4) || isempty(tol)
    if isa(A,'double')
        tol = 4e-292; % realmax('double')/eps('single')
    else 
        tol = 7e-31; % realmax('single')/eps('single')
    end
end
if (nargin < 3) || isempty(econ)
    econ = 1;
end
if size(J,1) == size(J,2)
    if size(J,1) ~= m
        error('J must have size(A,1)*size(A,1)')
    end
    J = sign2(diag(J)); % keep only diagonal entries
else
    if length(J) ~= m
        error('J must have size(A,1) entries')
    end
    J = sign2(J);
end
if any(J==0)
    error('J must be a signature matrix')
end
if issparse(A)
    error('JQR cannot handle sparse matrices')
end
od_method = true; % use orthogonal-diagonal method for real matrices


% initialize matrices
R = A;      % upper-triangular matrix
Q = eye(m); % J-orthogonal matrix
perm = 1:m; % row permutations

% start rotations
r = m+1; % row index
c = n+1; % column index

if isreal(R) % real Givens rotations
    while c > max(1,n-m+1);
        c = c-1;
        r = r-1;
        r1 = r;
        r2 = r1;
        retry = 0;
        r2_saved = [];
        while (r2 > 1 || retry > 0) % cancel entries in the column
            if r2 > 1 % proceed to cancel entries
                r2 = r2-1;
            else % retry to cancel saved entries
                r2 = r2_saved(retry);
                retry = retry - 1;
            end
            a1 = R(r1,c);
            a2 = R(r2,c);
            if abs(a2) <= tol*abs(a1)
                t = sign2(a1);
                R(r1,:) = t.*R(r1,:);
                R(r2,:) = t.*R(r2,:);
                R(r2,c) = 0;
                Q(r1,:) = t.*Q(r1,:);
                Q(r2,:) = t.*Q(r2,:);
                
            else % non-zero a2
                if J(r1) == J(r2) % same sign2, apply orthogonal rotation
                    
                    % apply Givens rotation using mixed-downdating (from [3])
                    if abs(a1) <= tol*abs(a2)
                        t = sign2(a2);
                        s = R(r1,:);
                        R(r1,:) = t.*R(r2,:);
                        R(r2,:) = -t.*s;
                        R(r2,c) = 0;
                        q = Q(r1,:);
                        Q(r1,:) = t.*Q(r2,:);
                        Q(r2,:) = -t.*q;
                        
                    elseif  abs(a1) > abs(a2)
                        t = a2/a1;
                        z = sign2(a1)*sqrt(1+t^2);
                        R(r1,:) = (R(r1,:) + t.*R(r2,:))./z;
                        R(r2,:) = -t.*R(r1,:) + z.*R(r2,:);
                        R(r2,c) = 0;
                        Q(r1,:) = (Q(r1,:) + t.*Q(r2,:))./z;
                        Q(r2,:) = -t.*Q(r1,:) + z.*Q(r2,:);
                        
                    else
                        t = a1/a2;
                        z = sign2(a2)*sqrt(1+t^2);
                        R(r1,:) = (t.*R(r1,:) + R(r2,:))./z;
                        R(r2,:) = (-R(r1,:) + z.*R(r2,:))./t;
                        R(r2,c) = 0;
                        Q(r1,:) = (t.*Q(r1,:) + Q(r2,:))./z;
                        Q(r2,:) = (-Q(r1,:) + z.*Q(r2,:))./t;
                    end
                    
                else % different sign2, apply hyperbolic rotation..
                    
                    da = abs(a1)-abs(a2);
                    if abs(da) > tol*(abs(a1)+abs(a2))
                        if da < 0 % permute rows if necessary
                            ind1 = [r1 r2];
                            ind2 = [r2 r1];
                            Q(ind2,:) = Q(ind1,:);
                            Q(:,ind2) = Q(:,ind1);
                            R(ind2,:) = R(ind1,:);
                            J(ind2) = J(ind1);
                            perm(ind2) = perm(ind1);
                            a3 = a1;
                            a1 = a2;
                            a2 = a3;
                        end
                        
                        if abs(a2) <= tol*abs(a1)
                            t = sign2(a1);
                            R(r1,:) = t.*R(r1,:);
                            R(r2,:) = t.*R(r2,:);
                            R(r2,c) = 0;
                            Q(r1,:) = t.*Q(r1,:);
                            Q(r2,:) = t.*Q(r2,:);
                        else % non-zero a2
                            if od_method
                                % apply hyperbolic rotation using orthogonal-diagonal
                                % method (forward numerical stable) (from [2])
                                a = 0.5*sign2(a1);
                                t = a.*sqrt((a1+a2)/(a1-a2));
                                z = a.*sqrt((a1-a2)/(a1+a2));
                                x = R(r1,:) - R(r2,:);
                                y = R(r1,:) + R(r2,:);
                                x = t.*x;
                                y = z.*y;
                                R(r1,:) = x + y;
                                R(r2,:) = y - x;
                                R(r2,c) = 0;
                                x = Q(r1,:) - Q(r2,:);
                                y = Q(r1,:) + Q(r2,:);
                                x = t.*x;
                                y = z.*y;
                                Q(r1,:) = x + y;
                                Q(r2,:) = y - x;
                            else
                                % apply hyperbolic rotation using mixed-downdating
                                t = a2/a1;
                                z = sign2(a1)*sqrt(1-t^2);
                                R(r1,:) = (R(r1,:) - t.*R(r2,:))./z;
                                R(r2,:) = -t.*R(r1,:) + z.*R(r2,:);
                                R(r2,c) = 0;
                                Q(r1,:) = (Q(r1,:) - t.*Q(r2,:))./z;
                                Q(r2,:) = -t.*Q(r1,:) + z.*Q(r2,:);
                            end
                        end
                        
                    else % no hyperbolic rotation
                        if (r2 > 1 || retry > 0) % still entries remaining
                            % cancelling all remaining entries first
                            retry = retry + 1;
                            r2_saved = [r2 r2_saved];
                            
                        else % no hyperbolic rotation
                            error('JQL failed to apply hyperbolic rotation. Check tolerance or method.')
                        end
                    end
                end
            end
        end
    end
    
    % make last diagonal entry positive
    if R(r,c) < 0
        Q(r,:) = -Q(r,:);
        R(r,:) = -R(r,:);
    end
    
else % complex Givens rotations
    while c > max(1,n-m+1);
        c = c-1;
        r = r-1;
        r1 = r;
        r2 = r1;
        retry = 0;
        r2_saved = [];
        while (r2 > 1 || retry > 0) % cancel entries in the column
            if r2 > 1 % proceed to cancel entries
                r2 = r2-1;
            else % retry to cancel saved entries
                r2 = r2_saved(retry);
                retry = retry - 1;
            end
            a1 = R(r1,c);
            a2 = R(r2,c);
            if abs(a2) <= tol*abs(a1)
                t = sign2(a1);
                R(r1,:) = t.*R(r1,:);
                R(r2,:) = t.*R(r2,:);
                R(r2,c) = 0;
                Q(r1,:) = t.*Q(r1,:);
                Q(r2,:) = t.*Q(r2,:);
                
            else % non-zero a2
                if J(r1) == J(r2) % same sign2, apply orthogonal rotation
                    
                    % apply Givens rotation using mixed-downdating (from [3])
                    if abs(a1) <= tol*abs(a2)
                        t = sign2(conj(a2));
                        s = R(r1,:);
                        R(r1,:) = t.*R(r2,:);
                        R(r2,:) = -conj(t).*s;
                        R(r2,c) = 0;
                        q = Q(r1,:);
                        Q(r1,:) = t.*Q(r2,:);
                        Q(r2,:) = -conj(t).*q;
                    else
                        a1a = abs(real(a1)) + abs(imag(a1));
                        a2a = abs(real(a2)) + abs(imag(a2));
                        
                        if  a1a > a2a
                            a1s = a1/a1a;
                            a11d = real(a1s)^2 + imag(a1s)^2;
                            a21s = a2/a1a;
                            a21d = real(a21s)^2 + imag(a21s)^2;
                            z = 1/(sign2(real(a1))*sqrt(1+a21d/a11d));
                            t = (conj(a21s)*a1s)/a11d;
                        else
                            a2s = a2/a2a;
                            a12s = a1/a2a;
                            a12d = real(a12s)^2 + imag(a12s)^2;
                            a22d = real(a2s)^2 + imag(a2s)^2;
                            z = abs(a1)/(sign2(real(a1))*a2a*sqrt(a12d + a22d));
                            t = sign2(a1)*conj(a2)/abs(a1);
                        end
                        
                        s = R(r1,:);
                        R(r1,:) = (R(r1,:) + t.*R(r2,:)).*z;
                        R(r2,:) = -conj(t*z).*s + R(r2,:).*z;
                        R(r2,c) = 0;
                        q = Q(r1,:);
                        Q(r1,:) = (Q(r1,:) + t.*Q(r2,:)).*z;
                        Q(r2,:) = -conj(t*z).*q + Q(r2,:).*z;
                    end
                    
                else % different sign2, apply hyperbolic rotation..
                    
                    da = abs(a1)-abs(a2);
                    if abs(da) > tol*(abs(a1)+abs(a2))
                        if da < 0 % permute rows if necessary
                            ind1 = [r1 r2];
                            ind2 = [r2 r1];
                            Q(ind2,:) = Q(ind1,:);
                            Q(:,ind2) = Q(:,ind1);
                            R(ind2,:) = R(ind1,:);
                            J(ind2) = J(ind1);
                            perm(ind2) = perm(ind1);
                            a3 = a1;
                            a1 = a2;
                            a2 = a3;
                        end
                        
                        if abs(a2) <= tol*abs(a1)
                            t = sign2(a1);
                            R(r1,:) = t.*R(r1,:);
                            R(r2,:) = t.*R(r2,:);
                            R(r2,c) = 0;
                            Q(r1,:) = t.*Q(r1,:);
                            Q(r2,:) = t.*Q(r2,:);
                        else % non-zero a2
                            % apply hyperbolic rotation using mixed-downdating
                            a1a = abs(real(a1)) + abs(imag(a1));
                            a2a = abs(real(a2)) + abs(imag(a2));
                            
                            if  a1a > a2a
                                a1s = a1/a1a;
                                a11d = real(a1s)^2 + imag(a1s)^2;
                                a21s = a2/a1a;
                                a21d = real(a21s)^2 + imag(a21s)^2;
                                z = 1/(sign2(real(a1))*sqrt(1-a21d/a11d));
                                t = (conj(a21s)*a1s)/a11d;
                            else
                                a2s = a2/a2a;
                                a12s = a1/a2a;
                                a12d = real(a12s)^2 + imag(a12s)^2;
                                a22d = real(a2s)^2 + imag(a2s)^2;
                                z = abs(a1)/(sign2(real(a1))*a2a*sqrt(a12d - a22d));
                                t = sign2(a1)*conj(a2)/abs(a1);
                            end
                            
                            s = R(r1,:);
                            R(r1,:) = (R(r1,:) - t.*R(r2,:)).*z;
                            R(r2,:) = -conj(t*z).*s + R(r2,:).*z;
                            R(r2,c) = 0;
                            q = Q(r1,:);
                            Q(r1,:) = (Q(r1,:) - t.*Q(r2,:)).*z;
                            Q(r2,:) = -conj(t*z).*q + Q(r2,:).*z;
                        end
                        
                    else % no hyperbolic rotation
                        if (r2 > 1 || retry > 0) % still entries remaining
                            % cancelling all remaining entries first
                            retry = retry + 1;
                            r2_saved = [r2 r2_saved];
                            
                        else % no hyperbolic rotation
                            error('JQL failed to apply hyperbolic rotation. Check tolerance.')
                        end
                    end
                end
            end
        end
    end
    
    % make last diagonal entry positive
    if real(R(r,c)) < 0
        Q(r,:) = -Q(r,:);
        R(r,:) = -R(r,:);
    end
end

% post-processing
Q(perm,:) = diag(J)*Q'*diag(J);
J = diag(J);
if ~econ && m > n
    R = R(end-n+1:end,1:n);
    Q = Q(:,end-n+1:end);
    J = J(end-n+1:end,end-n+1:end);
end
end

% sign with only positive or negative return (no zero)
function y = sign2(x)
if x == 0
    y=1;
else
    y=sign(x);
end
end

