function Sender=findSender(Sensors,Model,Receiver)
   
    Sender=[];
 
    n=Model.n;
 
    for i=1:n

        if (Sensors(i).MCH==Receiver & Sensors(i).id~=Receiver)
            Sender=[Sender,Sensors(i).id]; %#ok
        end

    end 

end