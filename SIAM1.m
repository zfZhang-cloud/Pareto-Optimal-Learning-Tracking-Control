clc;
clear;


t=1:1:50;
y1_d=sin(0.1*pi*t);
y2_d=cos(0.1*pi*t)-(1./(t+1));
y3_d=sin(0.1*pi*t)+1/2*cos(0.1*pi*t)-1/2;



%%
global G_1 G_2 G_3
syms t
A=[0.1-0.05*sin(0.2*t),0,0.01*t;0.1,-0.01*t,-0.02*cos(0.5*t);0,0.4,0.1*cos(0.2*t)];
B=[1-0.1*sin((pi*t)^2),0;0.1*sin(pi*t),0.1;0.2,0.5+0.02*cos((pi*t)^2)];
C=[-2,1.2,2;1,0.04,2;-0.4,0.3,3];
class_designA=eye(150,150);
for i=1:49
    for j=i+1:50

        class_designA((j-1)*3+1:j*3,(i-1)*3+1:i*3)=subs(A,t,j-1)*class_designA((j-2)*3+1:(j-1)*3,(i-1)*3+1:i*3); 
    end
end
G=zeros(150,100);
for i=2:50
    for j=1:i-1

        G((i-1)*3+1:i*3,(j-1)*2+1:j*2)=C*class_designA((i-1)*3+1:i*3,(j-1)*3+1:j*3)*subs(B,t,j-1); 
    end
end
for i=1:50
        G((i-1)*3+1:i*3,(i-1)*2+1:i*2)=C*subs(B,t,i-1);
end
G_1=G(1:3:148,:);
G_2=G(2:3:149,:);
G_3=G(3:3:150,:);




%% u1=0,equal weighted
u(1:100,1)=zeros(100,1);
for k=1: 11
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
    %d=a_1(:,c);
    %l=size(d,2);
    u(1:100,k+1)=u(1:100,k)+(1/k)*a_1*[1/3;1/3;1/3];
                y1(1:50,k+1)=G_1*u(1:100,k+1);
    y2(1:50,k+1)=G_2*u(1:100,k+1);
    y3(1:50,k+1)=G_3*u(1:100,k+1);

 
end
figure(1)
hold on
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(2)
hold on
plot(i,J2(1:10))
figure(3)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
 set(0,'defaultlinelinewidth',2)
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(4)
i=1:10;
hold on
plot(i,prec,'s:g')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w3(k)=f;
    
end
w123=w1+w2+w3;
figure(5)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(6)
i=1:10;
hold on
plot(i,w)


%% u1=0,unequal weighted
u(1:100,1)=zeros(100,1);
for k=1: 11
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
    %d=a_1(:,c);
    %l=size(d,2);
    u(1:100,k+1)=u(1:100,k)+(1/k)*a_1*[0.5;0.4;0.1];
                y1(1:50,k+1)=G_1*u(1:100,k+1);
    y2(1:50,k+1)=G_2*u(1:100,k+1);
    y3(1:50,k+1)=G_3*u(1:100,k+1);

 
end
figure(1)
hold on
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(2)
hold on
plot(i,J2(1:10))
figure(3)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
 set(0,'defaultlinelinewidth',2)
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(4)
i=1:10;
hold on
plot(i,prec,'s:g')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w3(k)=f;
    
end
w123=w1+w2+w3;
figure(5)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(6)
i=1:10;
hold on
plot(i,w)
%%  u1=0,update law21
u(1:100,1)=zeros(100,1);
sigma=zeros(1,11);
w=zeros(100,11);
for k=1: 31
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
      
        %l=size(d,2);
        %x0=unifrnd(0,1,[1 l])';
        %x0=x0/sum(x0);
        %A=ones(1,l);
        %b=1;
        %lb=zeros(l,1);
        %option=optimset; option.LargeScale='off';option.Display='off';
        %[x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        [xx,yy]=kkt1(a_1);
        w(:,k)=vpa(xx);
        sigma(k)=vpa(yy);
        %sigma(k)=min([((G_1*w(:,k))'*(y1_d'-y1(1:50,k)))/((G_1*w(:,k))'*(G_1*w(:,k))),((G_2*w(:,k))'*(y2_d'-y2(1:50,k)))/((G_2*w(:,k))'*(G_2*w(:,k))),((G_3*w(:,k))'*(y3_d'-y3(1:50,k)))/((G_3*w(:,k))'*(G_3*w(:,k)))]);
        u(1:100,k+1)=u(1:100,k)+sigma(k)*w(:,k);
end


figure(1)
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(2)
hold on
plot(i,J2(1:10))
figure(3)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(4)
i=1:10;
hold on
plot(i,prec,'s-')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
w123=w1+w2+w3;
figure(5)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(6)
i=1:10;
hold on
plot(i,w)














%%  %% u1=0,Algorithm4
u(1:100,1)=zeros(100,1);
w=zeros(100,30);
sigma=zeros(1,30);
for k=1: 11
    for i=1:3
        y1(1:50,3*(k-1)+i)=G_1*u(1:100,3*(k-1)+i);
        y2(1:50,3*(k-1)+i)=G_2*u(1:100,3*(k-1)+i);
        y3(1:50,3*(k-1)+i)=G_3*u(1:100,3*(k-1)+i);
        J1(3*(k-1)+i)=(y1_d'-y1(1:50,3*(k-1)+i))'*(y1_d'-y1(1:50,3*(k-1)+i));
        J2(3*(k-1)+i)=(y2_d'-y2(1:50,3*(k-1)+i))'*(y2_d'-y2(1:50,3*(k-1)+i));
        J3(3*(k-1)+i)=(y3_d'-y3(1:50,3*(k-1)+i))'*(y3_d'-y3(1:50,3*(k-1)+i));
        %F(3*(k-1)+i)=max([J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)]);
        %a=[J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)];
        a_1=[G_1'*(y1_d'-y1(1:50,3*(k-1)+i)),G_2'*(y2_d'-y2(1:50,3*(k-1)+i)),G_3'*(y3_d'-y3(1:50,3*(k-1)+i))];
        %b_1=sort([J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)]);
        %c=find(a>=b_1(4-i));
        P=[G_1',G_2',G_3']';
        P=P((i-1)*50+1:i*50,:);
        global d 
     a_1(:,i)=[];
        gamma=eye(100)-pinv(P)*P;
                d=gamma*a_1;
        %=size(d,2);
        %x0=unifrnd(0,1,[1 l])';
        %x0=x0/sum(x0);
        %A=ones(1,l);
        %b=1;
        %lb=zeros(l,1);
        %option=optimset; option.LargeScale='off';option.Display='off';
        %[x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        %global v1 e_1 e_2 e_3;
        %v1=d*x;
        if i==1
            w(:,k)=vpa(kkt12(d));
        end
        if i==2
             w(:,k)=vpa(kkt12(d));
        end
       if i==3
             w(:,k)=vpa(kkt12(d));
       end
        e_1=(y1_d'-y1(1:50,3*(k-1)+i));
        e_2=(y2_d'-y2(1:50,3*(k-1)+i));
        e_3=(y3_d'-y3(1:50,3*(k-1)+i));
        step=[((G_1*w(:,k))'*e_1)/((G_1*w(:,k))'*(G_1*w(:,k))),((G_2*w(:,k))'*e_2)/((G_2*w(:,k))'*(G_2*w(:,k))),((G_3*w(:,k))'*e_3)/((G_3*w(:,k))'*(G_3*w(:,k)))];
        step(:,i)=[];
        sigma(k)=min(step);
            u(1:100,3*(k-1)+i+1)=u(1:100,3*(k-1)+i)+sigma(k)*w(:,k);
           
    end
end

figure(1)
i=1:1:10;
hold on
plot(i,J1(1:3:30))
figure(2)
hold on
plot(i,J2(1:3:30))
figure(3)
hold on
plot(i,J3(1:3:30))
%»­Í¼ÊÇ·ñÖ§Åä
prec=zeros(1,10);
for k=1:10
    if J1(3*(k-1)+1)>=J1(3*(k)+1)&& J2(3*(k-1)+1)>=J2(3*(k)+1)&& J3(3*(k-1)+1)>=J3(3*(k)+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(4)
i=1:10;
hold on
plot(i,prec,'y<-')

w1=zeros(1,10);
for k=1:10
    a_1=[G_2'*(y2_d'-y2(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        v1(:,k)=vpa(kkt12(d));
        w1(k)=v1(:,k)'*v1(:,k);
    
end
%w2
w2=zeros(1,10);
for k=1:10
    a_2=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_2;
        v2(:,k)=vpa(kkt12(d));
        w2(k)=v2(:,k)'*v2(:,k);
    
end
%w3
w3=zeros(1,10);
for k=1:10
    a_3=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_2'*(y2_d'-y2(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_3;
        v3(:,k)=vpa(kkt12(d));
        w3(k)=v3(:,k)'*v3(:,k);
end
w123=w1+w2+w3;
figure(5)
i=1:10;
hold on
plot(i,w123)


w=zeros(1,10);
for k=1:10
     
    a_1=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_2'*(y2_d'-y2(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end
figure(6)
i=1:10;
hold on
plot(i,w)

%% u1=P_1{+}y_1d,equal weighted
u(1:100,1)=pinv(G_1)*y1_d';
for k=1: 11
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
    %d=a_1(:,c);
    %l=size(d,2);
    u(1:100,k+1)=u(1:100,k)+(1/k)*a_1*[1/3;1/3;1/3];
                y1(1:50,k+1)=G_1*u(1:100,k+1);
    y2(1:50,k+1)=G_2*u(1:100,k+1);
    y3(1:50,k+1)=G_3*u(1:100,k+1);

 
end
figure(11)
hold on
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(12)
hold on
plot(i,J2(1:10))
figure(13)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
 set(0,'defaultlinelinewidth',2)
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(14)
i=1:10;
hold on
plot(i,prec,'s:g')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w3(k)=f;
    
end
w123=w1+w2+w3;
figure(15)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(16)
i=1:10;
hold on
plot(i,w)


%% u1=P_1{+}y_1d,unequal weighted
u(1:100,1)=pinv(G_1)*y1_d';
for k=1: 11
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
    %d=a_1(:,c);
    %l=size(d,2);
    u(1:100,k+1)=u(1:100,k)+(1/k)*a_1*[0.5;0.4;0.1];
                y1(1:50,k+1)=G_1*u(1:100,k+1);
    y2(1:50,k+1)=G_2*u(1:100,k+1);
    y3(1:50,k+1)=G_3*u(1:100,k+1);

 
end
figure(11)
hold on
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(12)
hold on
plot(i,J2(1:10))
figure(13)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
 set(0,'defaultlinelinewidth',2)
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(14)
i=1:10;
hold on
plot(i,prec,'s:g')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w3(k)=f;
    
end
w123=w1+w2+w3;
figure(15)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(16)
i=1:10;
hold on
plot(i,w)
%%  u1=P_1{+}y_1d,update law21
u(1:100,1)=pinv(G_1)*y1_d';
sigma=zeros(1,11);
w=zeros(100,11);
for k=1: 31
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    %[F(k),index(k)]=max([J1(k),J2(k),J3(k)]);
    %F(k)=max([J1(k),J2(k),J3(k)]);
    %a=[J1(k),J2(k),J3(k)];
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
      
        %l=size(d,2);
        %x0=unifrnd(0,1,[1 l])';
        %x0=x0/sum(x0);
        %A=ones(1,l);
        %b=1;
        %lb=zeros(l,1);
        %option=optimset; option.LargeScale='off';option.Display='off';
        %[x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        [xx,yy]=kkt1(a_1);
        w(:,k)=vpa(xx);
        sigma(k)=vpa(yy);
        %sigma(k)=min([((G_1*w(:,k))'*(y1_d'-y1(1:50,k)))/((G_1*w(:,k))'*(G_1*w(:,k))),((G_2*w(:,k))'*(y2_d'-y2(1:50,k)))/((G_2*w(:,k))'*(G_2*w(:,k))),((G_3*w(:,k))'*(y3_d'-y3(1:50,k)))/((G_3*w(:,k))'*(G_3*w(:,k)))]);
        u(1:100,k+1)=u(1:100,k)+sigma(k)*w(:,k);
end


figure(11)
i=1:1:10;
hold on
plot(i,J1(1:10))
figure(12)
hold on
plot(i,J2(1:10))
figure(13)
hold on
plot(i,J3(1:10))
%»­Í¼ÊÇ·ñÖ§Åä
prec=zeros(1,10);
for k=1:10
    if J1(k)>=J1(k+1)&& J2(k)>=J2(k+1)&& J3(k)>=J3(k+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(14)
i=1:10;
hold on
plot(i,prec,'s-')


%¼ÆËãwi
w1=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w1(k)=f;
    
end
%w2
w2=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
%w3
w3=zeros(1,10);
for k=1:10
    y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w2(k)=f;
    
end
w123=w1+w2+w3;
figure(15)
i=1:10;
hold on
plot(i,w123)

%¼ÆËãw

w=zeros(1,10);
for k=1:10
     y1(1:50,k)=G_1*u(1:100,k);
    y2(1:50,k)=G_2*u(1:100,k);
    y3(1:50,k)=G_3*u(1:100,k);
    J1(k)=(y1_d'-y1(1:50,k))'*(y1_d'-y1(1:50,k));
    J2(k)=(y2_d'-y2(1:50,k))'*(y2_d'-y2(1:50,k));
    J3(k)=(y3_d'-y3(1:50,k))'*(y3_d'-y3(1:50,k));
    a_1=[G_1'*(y1_d'-y1(1:50,k)),G_2'*(y2_d'-y2(1:50,k)),G_3'*(y3_d'-y3(1:50,k))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end


figure(16)
i=1:10;
hold on
plot(i,w)














%%  %%u1=P_1{+}y_1d,Algorithm4
u(1:100,1)=pinv(G_1)*y1_d';
w=zeros(100,30);
sigma=zeros(1,30);
for k=1: 11
    for i=1:3
        y1(1:50,3*(k-1)+i)=G_1*u(1:100,3*(k-1)+i);
        y2(1:50,3*(k-1)+i)=G_2*u(1:100,3*(k-1)+i);
        y3(1:50,3*(k-1)+i)=G_3*u(1:100,3*(k-1)+i);
        J1(3*(k-1)+i)=(y1_d'-y1(1:50,3*(k-1)+i))'*(y1_d'-y1(1:50,3*(k-1)+i));
        J2(3*(k-1)+i)=(y2_d'-y2(1:50,3*(k-1)+i))'*(y2_d'-y2(1:50,3*(k-1)+i));
        J3(3*(k-1)+i)=(y3_d'-y3(1:50,3*(k-1)+i))'*(y3_d'-y3(1:50,3*(k-1)+i));
        %F(3*(k-1)+i)=max([J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)]);
        %a=[J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)];
        a_1=[G_1'*(y1_d'-y1(1:50,3*(k-1)+i)),G_2'*(y2_d'-y2(1:50,3*(k-1)+i)),G_3'*(y3_d'-y3(1:50,3*(k-1)+i))];
        %b_1=sort([J1(3*(k-1)+i),J2(3*(k-1)+i),J3(3*(k-1)+i)]);
        %c=find(a>=b_1(4-i));
        P=[G_1',G_2',G_3']';
        P=P((i-1)*50+1:i*50,:);
        global d 
     a_1(:,i)=[];
        gamma=eye(100)-pinv(P)*P;
                d=gamma*a_1;
        %=size(d,2);
        %x0=unifrnd(0,1,[1 l])';
        %x0=x0/sum(x0);
        %A=ones(1,l);
        %b=1;
        %lb=zeros(l,1);
        %option=optimset; option.LargeScale='off';option.Display='off';
        %[x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        %global v1 e_1 e_2 e_3;
        %v1=d*x;
        if i==1
            w(:,k)=vpa(kkt12(d));
        end
        if i==2
             w(:,k)=vpa(kkt12(d));
        end
       if i==3
             w(:,k)=vpa(kkt12(d));
       end
        e_1=(y1_d'-y1(1:50,3*(k-1)+i));
        e_2=(y2_d'-y2(1:50,3*(k-1)+i));
        e_3=(y3_d'-y3(1:50,3*(k-1)+i));
        step=[((G_1*w(:,k))'*e_1)/((G_1*w(:,k))'*(G_1*w(:,k))),((G_2*w(:,k))'*e_2)/((G_2*w(:,k))'*(G_2*w(:,k))),((G_3*w(:,k))'*e_3)/((G_3*w(:,k))'*(G_3*w(:,k)))];
        step(:,i)=[];
        sigma(k)=min(step);
            u(1:100,3*(k-1)+i+1)=u(1:100,3*(k-1)+i)+sigma(k)*w(:,k);
           
    end
end

figure(11)
i=1:1:10;
hold on
plot(i,J1(1:3:30))
figure(12)
hold on
plot(i,J2(1:3:30))
figure(13)
hold on
plot(i,J3(1:3:30))
%»­Í¼ÊÇ·ñÖ§Åä
prec=zeros(1,10);
for k=1:10
    if J1(3*(k-1)+1)>=J1(3*(k)+1)&& J2(3*(k-1)+1)>=J2(3*(k)+1)&& J3(3*(k-1)+1)>=J3(3*(k)+1)
        if J1(k)==J1(k+1) && J2(k)==J2(k+1)&& J3(k)==J3(k+1)
            prec(k)=1/2;
        else
            prec(k)=1;
        end
    else
        prec(k)=0;
    end
end
figure(14)
i=1:10;
hold on
plot(i,prec,'y<-')

w1=zeros(1,10);
for k=1:10
    a_1=[G_2'*(y2_d'-y2(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_1)*G_1;
                d=gamma*a_1;
        v1(:,k)=vpa(kkt12(d));
        w1(k)=v1(:,k)'*v1(:,k);
    
end
%w2
w2=zeros(1,10);
for k=1:10
    a_2=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_2)*G_2;
                d=gamma*a_2;
        v2(:,k)=vpa(kkt12(d));
        w2(k)=v2(:,k)'*v2(:,k);
    
end
%w3
w3=zeros(1,10);
for k=1:10
    a_3=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_2'*(y2_d'-y2(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        gamma=eye(100)-pinv(G_3)*G_3;
                d=gamma*a_3;
        v3(:,k)=vpa(kkt12(d));
        w3(k)=v3(:,k)'*v3(:,k);
end
w123=w1+w2+w3;
figure(15)
i=1:10;
hold on
plot(i,w123)


w=zeros(1,10);
for k=1:10
     
    a_1=[G_1'*(y1_d'-y1(1:50,3*(k-1)+1)),G_2'*(y2_d'-y2(1:50,3*(k-1)+1)),G_3'*(y3_d'-y3(1:50,3*(k-1)+1))];
    %b_1=sort([J1(k),J2(k),J3(k)]);
    %c=find(a>=b_1(3));
        global d 
        d=a_1;
        l=size(d,2);
        x0=unifrnd(0,1,[1 l])';
        x0=x0/sum(x0);
        A=ones(1,l);
        b=1;
        lb=zeros(l,1);
        option=optimset; option.LargeScale='off';option.Display='off';
        [x,f]=fmincon('my_sum',x0,[],[],A,b,lb,[],[],option);
        w(k)=f;
    
end
figure(16)
i=1:10;
hold on
plot(i,w)

