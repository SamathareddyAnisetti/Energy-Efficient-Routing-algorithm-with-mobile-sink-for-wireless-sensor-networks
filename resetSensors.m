function Sensors=resetSensors(Sensors,Model)

    n=Model.n;
    for i=1:n
        Sensors(i).MCH=n+1;
        Sensors(i).type='N';
        Sensors(i).dis2ch=inf;
    end
    
end