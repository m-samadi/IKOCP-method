function RWP_Point = RWP_Point_Set(RWP_Point,RWP_Point_Count,Network_Length,Network_Width)

    for i=1:RWP_Point_Count
        RWP_Point(i,2)=ceil(rand()*Network_Length);
        RWP_Point(i,3)=ceil(rand()*Network_Width);
    end;
    
end

