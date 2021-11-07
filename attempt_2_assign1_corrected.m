clear all
clc
x_nodes = 11;
y_nodes = 11;
x=linspace(0,1,x_nodes); %Grid
y=linspace(0,1,y_nodes);
q = 1353;
k = 0.69;
h = 8;
L = 0.5;
l = 0.1;
d_x=L/(x_nodes-1); %Distance between two grid points in x and y
d_y=L/(y_nodes-1);
beta = d_x/d_y;
T = zeros(x_nodes,y_nodes);
no_of_iter = 0;
for i=4:x_nodes
    for j=4:y_nodes
        T(i,j)=30;
    end
end
Told = T;
error=9e9;
%Setting the tolerance value (This affects the accuracy and time for
%convergence)
tol=1e-06;

while(error > tol)
    %interior nodes
    for i=2
        for j=2:y_nodes-1
            T(i,j) = (1/(2+2*beta*beta))*(Told(i-1,j) + Told(i+1,j) + ((beta*beta)*(Told(i,j-1)+Told(i,j+1))) + (d_x *d_x)*(q/k));
        end
    end
    for j=2
        for i=3:x_nodes-1
            T(i,j) = (1/(2+2*beta*beta))*(Told(i-1,j) + Told(i+1,j) + ((beta*beta)*(Told(i,j-1)+Told(i,j+1))) + (d_x *d_x)*(q/k));
        end
    end
    for i=3
        for j=3
            T(i,j) = (1/(2+2*beta*beta))*(Told(i-1,j) + Told(i+1,j) + ((beta*beta)*(Told(i,j-1)+Told(i,j+1))) + (d_x *d_x)*(q/k));
        end
    end
    %boundary 1-2
    for i=1
        for j=2:y_nodes-1
            T(i,j) = 0.25*(2*Told(i+1,j)+Told(i,j+1)+Told(i,j-1)+(d_x*d_x*q/k));
        end
    end
    
    %boundary 2-3
    for i=2
        for j=11
            T(i,j) = (1/(4+(2*h*d_x/k)))*(2*Told(i,j-1)+Told(i-1,j)+Told(i+1,j)+(60*h*d_x/k)+(d_x*d_x*q/k));
        end
    end
    %boundary 3-a
    for i=3
        for j=4:y_nodes-1
            T(i,j) = (1/(4+(2*h*d_x/k)))*(2*Told(i-1,j)+Told(i,j+1)+Told(i,j-1)+(60*h*d_x/k)+(d_x*d_x*q/k));
        end
    end
    %boundary a-4
    for j=3
        for i=4:x_nodes-1
            T(i,j) = (1/(4+(2*h*d_x/k)))*(2*Told(i,j-1)+Told(i+1,j)+Told(i-1,j)+(60*h*d_x/k)+(d_x*d_x*q/k));
        end
    end
    %boundary 4-5
    for i=11 
        for j=2
            T(i,j)= (1/(4+(2*h*d_x/k)))*(2*Told(i-1,j)+Told(i,j-1)+Told(i,j+1)+(60*h*d_x/k)+(d_x*d_x*q/k));
        end
    end
    %boundary 1-5
    for i=2:x_nodes-1
        for j=1
            T(i,j) =  0.25*(2*Told(i,j+1)+Told(i+1,j)+Told(i-1,j)+(d_x*d_x*q/k));
        end
    end
    %Corner point 1
    T(1,1)=0.25*(2*Told(2,1)+ 2*Told(1,2) + (d_x*d_x*q/k));
    %Corner point 2
    T(1,11)=(1/(4+(2*h*d_x/k)))*(2*Told(2,11)+2*Told(1,10)+(60*h*d_x/k)+(d_x*d_x*q/k));
    %Corner point 3
    T(3,11)=(1/(4+(4*h*d_x/k)))*(2*Told(2,11)+2*Told(3,10)+(120*h*d_x/k)+(d_x*d_x*q/k));
    %Corner point 4
    T(11,3)=(1/(4+(4*h*d_x/k)))*(2*Told(11,2)+2*Told(10,3)+(120*h*d_x/k)+(d_x*d_x*q/k));
    %Corner point 5
    T(11,1)=(1/(4+(2*h*d_x/k)))*(2*Told(11,2)+2*Told(10,1)+(60*h*d_x/k)+(d_x*d_x*q/k));
    
    error1=max(abs(Told-T)); 
    error=max(error1);%Obtains the final value of the error.
    Told=T; %Fixing the Told values to the latest iteration
    no_of_iter=no_of_iter+1%Incrementing the number of iterations
end
contourf(x,y,T,10,'showText','on');colormap(jet);
colorbar;
title(sprintf('Temperature Distribution at the end'));