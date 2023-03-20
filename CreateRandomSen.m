function CreateRandomSen(Model,Area)
    n=Model.n;
    x=Area.x;
    y=Area.y;
    X=zeros(1,n);
    Y=X;
    for i=1:1:n
        X(i)=rand()*x;
        Y(i)=rand()*y;
    end
    save ('Locations','X','Y');

end