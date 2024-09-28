function [rlist,plist] = add_atoms(xlist,ylist,xa,xb,Rb,N,Natoms)
    
    plist=zeros(Natoms,3);
    rlist=zeros(Natoms,1);
    count=0;
    while(count<Natoms)
        index=fix(rand()*N)+1; % choose the interval that we sample from 
        xtest=1-rand()*(1-xlist(index));
        ytest=ylist(index)+rand()*(ylist(index+1)-ylist(index));
        % count=count+1;
        % rlist(count)=xtest;
        if xtest>xlist(index+1)
            count=count+1;
            rlist(count)=xtest;
        elseif ytest>xa*xtest^4+xb*xtest^2
            count=count+1;
            rlist(count)=xtest;
        end
    end

    % generate the atoms coordinates for each shell using Archimedes' Theorem
    for i=1:Natoms
        theta=2*pi*rand();
        z=2*rand()-1;
        x=sqrt(1-z^2)*cos(theta);
        y=sqrt(1-z^2)*sin(theta);
        plist(i,:)=rlist(i)*[x,y,z];
    end
end