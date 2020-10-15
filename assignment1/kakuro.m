% Ϊ����ƣ�����ʹ���������� a �� b ���ֱ��ź����Լ�������������Լ��������
function [out,k]=kakuro(a,b)  
    if nargin<2 %û���ر������ʱ��ֱ��ѡ��������
        a=[inf inf inf inf inf inf inf inf inf inf;
           inf 3   0   0   4   0   0   inf inf inf;
           21  0   0   0   0   0   0   5   0   0  ;
           6   0   0   29  0   0   0   0   0   0  ;
           inf 26  0   0   0   0   13  0   0   0  ;
           15  0   0   0   11  0   0   0   0   inf;
           10  0   0   0   0   17  0   0   0   inf;        
           inf inf 15  0   0   inf 7   0   0   inf;];     
        b=[inf inf 21  6   inf 22  7   inf inf inf;        
           inf 3   0   0   21  0   0   inf 38  6  ;        
           inf 0   0   0   0   0   0   16  0   0  ;        
           inf 0   0   23  0   0   0   0   0   0  ;        
           inf 3   0   0   0   0   3   0   0   0  ;        
           inf 0   0   0   13  0   0   0   0   inf;        
           inf 0   0   0   0   inf 0   0   0   inf;        
           inf inf inf 0   0   inf inf 0   0   inf]; 
    end
    
    n=length(a(1,:)); 
    m=length(b(:,1)); 

    for t=1:n 
        c{t}=GenFromVec(b(:,t));%�г������еĿ����� 
    end
    for t=1:m     
        d{t}=GenFromVec(a(t,:));%�г������еĿ�����    
    end
    
    d=prune(d,c,n,m);
    c=prune(c,d,m,n);%����ɸѡ 
    
    Ans=solveProblem(c,d,1,{[]});%�ݹ���� 
    if length(Ans)==1               %������
        fprintf('No Answer') 
    else
        for i=1:(length(Ans)-1)         
            out(:,:,i)=Ans{i+1};     
        end
    end
end

function y=prune(d,c,n,m) %��֦��
    for h=1:m 
        for index=1:length(d{h}(:,1))         
            tmp=d{h}(index,:);                  
            for i=1:n            
                if  Check(tmp(i),c{i})~=1        %û�п��������һ����    
                    d{h}(index,:)=zeros(1,n);
                    break            
                end
            end
        end
        j=1;z=[]; 
        for index=1:length(d{h}(:,1))               %��ȥ��
            if d{h}(index,:)~=zeros(1,n)
                z(j,:)=d{h}(index,:);
                j=j+1;
            end
        end
        y{h}=z;
    end
end

function f=solveProblem(d,c,h,f) % ��������⡣
    n=length(c(1,:));
    m=length(d(1,:));
    e=c;
    for i1=1:length(d{h}(:,1))
        tmp=d{h}(i1,:);
        flag=1;
        c=e;
        for i=1:n           %����Ƿ���С�
            if isempty(c{i})~=1   
                 if  Check(tmp(i),c{i})~=1
                    flag=0;
                    break
                 end 
            else
                flag=0;
                break
            end
        end
        if flag==1
            for index1=1:n
                for index2=1:length(c{index1}(:,1))
                    if c{index1}(index2,h)~=tmp(index1)
                        c{index1}(index2,1)=0;
                    end
                end
                j=1;
                z=[];
                for index2=1:length(c{index1}(:,1))
                     if c{index1}(index2,1)~=0
                        z(j,:)=c{index1}(index2,:);
                        j=j+1;
                     end
                end
                c{index1}=z;
            end
            if h==m
                flag=1;
                for i=1:n
                    if isempty(c{i}) 
                        flag=0;
                        break
                    end
                end
                if flag==1
                    f_length=length(f);
                    for i=1:n
                        f{f_length+1}(i,:)=c{i};
                    end
                end 
            end
            if h<m  % ������һ��
                f= solveProblem(d,c,h+1,f);
            else
            end
        else
        end
    end
end

function z=GenFromVec(a)
    n=length(a);y=0;
    for i=1:n
        if a(i)~=inf&&a(i)~=0  % ��ʱ����һ����ʼΪ�����ո�
            isContinuous=1;
            p=0;j=i+1;
            while j<=n & a(j) == 0 
                p=p+1;
                j=j+1;
            end
            Divisions=GenMatrix(GenAll(a(i),p)); % �������п��еķֻ���
            y=VecCat(y,Divisions);
        end
    end
    
    z=zeros(1,n);
    for i=1:n
        if a(i)~=0
            z(1,i)=inf;
        end
    end
    
    m=length(y(:,1));   % Ϊÿ�ַֻ�����һ��������
    for i=2:m
        z(i,:)=z(1,:);
    end
    
    j=1;
    for i=1:length(z(1,:))
        if z(1,i)==0
            z(:,i)=y(:,j);
        j=j+1;
        end
    end
end

function y=Check(a,b) %�ж���ѡ�еĵ�n��Ԫ���Ƿ������n�еĿ�����
     n=length(b(1,:));
     m=length(b(:,1));
     y=0;
     for i=1:n
         for j=1:m
             if a==b(j,i)
                 y=1;
                 break
             end
         end
     end
end

function a=GenAll(b,c) %�������п��ܵķֽⲢ���ء�
    d=1:9;
    VecChosen=zeros(9,1);
    cnt=1;
    for i=0:1
        VecChosen(1)=i;
        for i=0:1
            VecChosen(2)=i;
            for i=0:1
                VecChosen(3)=i;
                for i=0:1
                    VecChosen(4)=i;
                    for i=0:1
                        VecChosen(5)=i;
                        for i=0:1
                            VecChosen(6)=i;
                            for i=0:1
                                VecChosen(7)=i;
                                for i=0:1
                                    VecChosen(8)=i;
                                    for i=0:1
                                        VecChosen(9)=i;
                                        if sum(VecChosen)==c&&d*VecChosen==b
                                            g=d.*(VecChosen)';
                                            a(cnt,:)=g(g~=0);
                                            cnt=cnt+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function a=GenMatrix(b)
    m=length(b(:,1));
    p=length(perms(b(1,:)));
    for i=1:m
        a(((i-1)*p+1:(i*p)),:)=perms(b(i,:));
    end
end

function c=VecCat(a,b) %������ strcat 
    if sum(a)==0
        c=b;
    else
        length_a=length(a(:,1));
        length_b=length(b(:,1));
        for i=1:length_a
            for j=1:length_b
                c((i-1)*length_b+j,:)=[a(i,:) b(j,:)]; %#ok<*AGROW>
            end
        end
    end
end