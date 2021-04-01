clc

A = randn(100);
J = eye(100);
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(100);
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(100,50);
J = eye(50);
[L,Q,Jp] = jlq(A,J,[],0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(100,50);
J = blkdiag(eye(20),-eye(30));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(50,100);
J = eye(100);
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(50,100);
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(100,50);
A = A*A';
J = eye(100);
[L,Q,Jp] = jlq(A,J);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = randn(100,50);
A = A*A';
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2



A = complex(randn(100),randn(100));
J = eye(100);
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(100),randn(100));
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(100,50),randn(100,50));
J = eye(50);
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(100,50),randn(100,50));
J = blkdiag(eye(20),-eye(30));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(50,100),randn(50,100));
J = eye(100);
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(50,100),randn(50,100));
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(100,50),randn(100,50));
A = A*A';
J = eye(100);
[L,Q,Jp] = jlq(A,J);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2

A = complex(randn(100,50),randn(100,50));
A = A*A';
J = blkdiag(eye(50),-eye(50));
[L,Q,Jp] = jlq(A,J,0);
norm(A-L*Q)
norm(Jp - Q*J*Q') > 1e-2