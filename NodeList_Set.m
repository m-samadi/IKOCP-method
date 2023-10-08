function NodeList = NodeList_Set(NodeList,Network_Length,Network_Width)

    [m n]=size(NodeList);
    for i=1:m
        NodeList(i,2)=ceil(rand()*Network_Length);  
        NodeList(i,3)=ceil(rand()*Network_Width);
    end;

end

