function P = lqr_p_list(A,B,Q,R,Qf,N)

[m,n] = size(Qf);
P = zeros(m,n,N);
P(:,:,N) = Qf;

for i=N-1:-1:1
    Pp = P(:,:,i+1);
    P(:,:,i) = Q + A'*Pp*A - A'*Pp*B/(R+B'*Pp*B)*B'*Pp*A;
end