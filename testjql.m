clc

A = randn(10);
J = eye(10);
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(10);
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(10,5);
J = eye(10);
[Q,L,Jp] = jql(A,J);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(10,5);
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(5,10);
J = eye(5);
[Q,L,Jp] = jql(A,J);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(5,10);
J = blkdiag(eye(2),-eye(3));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(10,5);
A = A*A';
J = eye(10);
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = randn(10,5);
A = A*A';
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2



A = complex(randn(10),randn(10));
J = eye(10);
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10),randn(10));
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
J = eye(10);
[Q,L,Jp] = jql(A,J);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(5,10),randn(5,10));
J = eye(5);
[Q,L,Jp] = jql(A,J);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(5,10),randn(5,10));
J = blkdiag(eye(2),-eye(3));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
A = A*A';
J = eye(10);
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2

A = complex(randn(10,5),randn(10,5));
A = A*A';
J = blkdiag(eye(5),-eye(5));
[Q,L,Jp] = jql(A,J,0);
norm(A-Q*L)
norm(Jp - Q'*J*Q) > 1e-2