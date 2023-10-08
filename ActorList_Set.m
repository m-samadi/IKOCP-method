function ActorList = ActorList_Set(ActorList,ActorNumber,RWP_Point_Count)

    for i=1:ActorNumber
        ActorList(i,6)=ceil(rand()*RWP_Point_Count); 
    end;

end

