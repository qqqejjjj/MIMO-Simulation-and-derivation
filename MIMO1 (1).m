rc=5;a=250;sum_of_inter=0;CDFof7=[];
for j=1:100
nodes=drawcelluar(rc,a);%nodes包含基站与待定用户信息
for i=1:length(nodes)
   r= sqrt(nodes(i).usersloc(1)^2+ nodes(i).usersloc(2)^2);
   rcenter2=nodes(i).centerloc(1)^2+ nodes(i).centerloc(2)^2;
    if ((r<=8*2*rc/sqrt(3))&&((abs(rcenter2-7*((2*rc)^2))<0.001)||(abs(rcenter2-7*4*((2*rc)^2))<0.001)||(abs(rcenter2-7*9*((2*rc)^2))<0.001)))||i==1
        z=10^(nodes(i).guass/10);
        r_r=r^3.8;
        nodes(i).betas=z/r_r;
        if i~=1
        sum_of_inter=sum_of_inter+(nodes(i).betas)^2;
        end
    end
end
SIRS=(nodes(1).betas)^2/sum_of_inter;
SIRS=10*log10(SIRS);
CDFof7(j)=SIRS;
end

rc=5;a=250;sum_of_inter=0;CDFof1=[];
for j=1:100
nodes=drawcelluar(rc,a);%nodes包含基站与待定用户信息
for i=1:length(nodes)
   r= sqrt(nodes(i).usersloc(1)^2+ nodes(i).usersloc(2)^2);
    if (r<=8*2*rc/sqrt(3))
        z=10^(nodes(i).guass/10);
        r_r=r^3.8;
        nodes(i).betas=z/r_r;
        if i~=1
        sum_of_inter=sum_of_inter+(nodes(i).betas)^2;
        end
    end
end

SIRS=(nodes(1).betas)^2/sum_of_inter;
SIRS=10*log10(SIRS);
CDFof1(j)=SIRS;
end

rc=5;a=250;sum_of_inter=0;CDFof3=[];
for j=1:100
nodes=drawcelluar(rc,a);%nodes包含基站与待定用户信息
for i=1:length(nodes)
   r= sqrt(nodes(i).usersloc(1)^2+ nodes(i).usersloc(2)^2);
   rcenter2=nodes(i).centerloc(1)^2+ nodes(i).centerloc(2)^2;
    if ((r<=8*2*rc/sqrt(3))&&((abs(rcenter2-3*((2*rc)^2))<0.001)||(abs(rcenter2-3*4*((2*rc)^2))<0.001)||(abs(rcenter2-3*9*((2*rc)^2))<0.001)))||i==1
        z=10^(nodes(i).guass/10);
        r_r=r^3.8;
        nodes(i).betas=z/r_r;
        if i~=1
        sum_of_inter=sum_of_inter+(nodes(i).betas)^2;
        end
    end
end
SIRS=(nodes(1).betas)^2/sum_of_inter;
SIRS=10*log10(SIRS);
CDFof3(j)=SIRS;
end

cdfplot(CDFof1);
hold on;
cdfplot(CDFof3);
hold on;
cdfplot(CDFof7);
xlabel('SIRS');
legend('1','3','7')


function nodes=drawcelluar(rc,a)

dy=2*rc;dx=rc*sqrt(3);
A = 1:7;
dots=[0,0];
A=pi/3*A;%圆周6等分
num=1;
for yk=[0:dy:a,0:-dy:-a]
    yfun=inline(['sqrt(3)*x/3+',num2str(yk)]);%内联函数
    for xk=[0:dx:a,0:-dx:-a]
        xp=xk;
        yp=yfun(xp);
        if -a/2<xp && xp<a/2 && -a/2<yp && yp<a/2
            T = (xp+1i*yp)+rc*exp(1i*A)*2/sqrt(3);%复数形式
%             Vertx = real(T)'+50;
%             Verty = imag(T)'+50;
%             V = [Vertx,Verty];%六个角点的坐标
           
            %fprintf('第%d个六边形加粗\n',num);
            %plot(T)
            %plot(0,0,'*b')
            hold on;
            %fprintf('%d %d\n',xp,yp);
            dot=cat(2,xp',yp');
            [Lia, ] = ismember(dot,dots, 'rows');
            if Lia==0||num==1
            dot1=[dots;dot];
            dots=dot1;
            num=num+1;
            i = 0;
            while i < 1
            x = 2*rc*rand(1,2)-1*rc;
            if (abs(x(1)) + abs(x(2))/sqrt(3) ) <= rc&& abs(x(2)) <= rc*sqrt(3)/2&& (abs(x(1)) + abs(x(2))/sqrt(3))  >= (rc/16)&& abs(x(2)) >= (rc/16)*sqrt(3)/2
            i = i+1;
            nodes(num-1).centerloc=dot;
            nodes(num-1).usersloc=[x(1)+xp,x(2)+yp];
            nodes(num-1).guass=normrnd(0,8)%生成方差为8的高斯分布  
            %plot(x(1) + xp, x(2) + yp,'r*'); 
            end
            end
            end
        end
    end
end
hold off
nodes;
end