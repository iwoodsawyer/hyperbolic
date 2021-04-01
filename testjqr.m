clc

A = rand(10);
J = eye(10);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(10);
J = blkdiag(-eye(5),eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(10,5);
J = eye(10);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(10,5);
J = blkdiag(eye(5),-eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(5,10);
J = eye(5);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(5,10);
J = blkdiag(eye(2),-eye(3));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(10,5);
A = A*A';
J = eye(10);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = rand(10,5);
A = A*A';
J = blkdiag(eye(5),-eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2


A = complex(randn(10),randn(10));
J = eye(10);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10),randn(10));
J = blkdiag(-eye(5),eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
J = eye(10);
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
J = blkdiag(eye(5),-eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(5,10),randn(5,10));
J = blkdiag(-eye(2),eye(3));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
if norm(Jp - Q'*J*Q) > 1e-2
 error('stop') 
end

A = complex(randn(5,10),randn(5,10));
J = blkdiag(eye(2),-eye(3));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
if norm(Jp - Q'*J*Q) > 1e-2
 error('stop') 
end,

A = complex(randn(10,5),randn(10,5));
A = A*A';
J = blkdiag(-eye(5),eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
if norm(Jp - Q'*J*Q) > 1e-2
 error('stop') 
end

A = complex(randn(10,5),randn(10,5));
A = A*A';
J = blkdiag(eye(5),-eye(5));
[Q,R,Jp] = jqr(A,J,0);
norm(A-Q*R)
if norm(Jp - Q'*J*Q) > 1e-2
 error('stop') 
end






















