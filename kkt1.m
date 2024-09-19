function [xx,yy]=kkt1(x)
global G_1 G_2 G_3
F=x'*x;
%case1
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=beta1;
equ5=beta2;
equ6=beta3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if lambda1>0&&lambda2>0&&lambda3>0
    xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case2
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=lambda1;
equ5=beta2;
equ6=beta3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if beta1>0&&lambda2>0&&lambda3>0
   xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case3
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=beta1;
equ5=lambda2;
equ6=beta3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if lambda1>0&&beta2>0&&lambda3>0
   xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case4
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=beta1;
equ5=beta2;
equ6=lambda3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if lambda1>0&&lambda2>0&&beta3>0
    xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case5
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=lambda1;
equ5=lambda2;
equ6=beta3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if beta1>0&&beta2>0&&lambda3>0
    xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case6
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=lambda1;
equ5=beta2;
equ6=lambda3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if beta1>0&&lambda2>0&&beta3>0
    xx=x*[lambda1;lambda2;lambda3];
   yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
%case7
syms lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1
equ1=2*F(1,:)*[lambda1;lambda2;lambda3]-alhpa1-beta1;
equ2=2*F(2,:)*[lambda1;lambda2;lambda3]-alhpa1-beta2;
equ3=2*F(3,:)*[lambda1;lambda2;lambda3]-alhpa1-beta3;
equ4=beta1;
equ5=lambda2;
equ6=lambda3;
equ7=lambda1+lambda2+lambda3-1;
[lambda1,lambda2,lambda3,beta1,beta2,beta3,alhpa1] = solve([equ1 equ2 equ3 equ4 equ5 equ6 equ7], [lambda1 lambda2 lambda3 beta1 beta2 beta3 alhpa1]);
if lambda1>0&&beta2>0&&beta3>0
    xx=x*[lambda1;lambda2;lambda3];
    yy=min([((x*[lambda1;lambda2;lambda3])'*x(:,1))/((G_1*(x*[lambda1;lambda2;lambda3]))'*(G_1*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,2))/((G_2*(x*[lambda1;lambda2;lambda3]))'*(G_2*(x*[lambda1;lambda2;lambda3]))),((x*[lambda1;lambda2;lambda3])'*x(:,3))/((G_3*(x*[lambda1;lambda2;lambda3]))'*(G_3*(x*[lambda1;lambda2;lambda3])))]);
    return
end
end