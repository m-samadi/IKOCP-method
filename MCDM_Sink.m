function Region = MCDM_Sink(Network_Length,Network_Width,NumberOfNodes)
    %%% Initialization
    SurfaceArea=Network_Length*Network_Width; 

    %%% Selected state
    % SurfaceArea   
    if (SurfaceArea>=0)&&(SurfaceArea<250000)
        Selected_SurfaceArea='Small';
    elseif (SurfaceArea>=250000)&&(SurfaceArea<1000000)    
        Selected_SurfaceArea='Mediocre';
    elseif (SurfaceArea>=1000000)&&(SurfaceArea<2250000)    
        Selected_SurfaceArea='Large'; 
    else
        Selected_SurfaceArea='VeryLarge';
    end;

    % NumberOfNodes    
    if (NumberOfNodes>=0)&&(NumberOfNodes<250)
        Selected_NumberOfNodes='Small';
    elseif (NumberOfNodes>=250)&&(NumberOfNodes<500)    
        Selected_NumberOfNodes='Medium';
    elseif (NumberOfNodes>=500)&&(NumberOfNodes<750)    
        Selected_NumberOfNodes='Large'; 
    else
        Selected_NumberOfNodes='VeryLarge';
    end; 

    % NumberOfRegions
    if (strcmp(Selected_SurfaceArea,'Small')==1)&&(strcmp(Selected_NumberOfNodes,'Small')==1)
        Selected_NumberOfRegions='Feeble';
    elseif (strcmp(Selected_SurfaceArea,'Small')==1)&&(strcmp(Selected_NumberOfNodes,'Medium')==1)
        Selected_NumberOfRegions='Feeble'; 
    elseif (strcmp(Selected_SurfaceArea,'Small')==1)&&(strcmp(Selected_NumberOfNodes,'Large')==1)
        Selected_NumberOfRegions='Feeble'; 
    elseif (strcmp(Selected_SurfaceArea,'Small')==1)&&(strcmp(Selected_NumberOfNodes,'VeryLarge')==1)
        Selected_NumberOfRegions='Few';     
    elseif (strcmp(Selected_SurfaceArea,'Mediocre')==1)&&(strcmp(Selected_NumberOfNodes,'Small')==1)
        Selected_NumberOfRegions='Few';
    elseif (strcmp(Selected_SurfaceArea,'Mediocre')==1)&&(strcmp(Selected_NumberOfNodes,'Medium')==1)
        Selected_NumberOfRegions='Few'; 
    elseif (strcmp(Selected_SurfaceArea,'Mediocre')==1)&&(strcmp(Selected_NumberOfNodes,'Large')==1)
        Selected_NumberOfRegions='Medium'; 
    elseif (strcmp(Selected_SurfaceArea,'Mediocre')==1)&&(strcmp(Selected_NumberOfNodes,'VeryLarge')==1)
        Selected_NumberOfRegions='Medium';      
    elseif (strcmp(Selected_SurfaceArea,'Large')==1)&&(strcmp(Selected_NumberOfNodes,'Small')==1)
        Selected_NumberOfRegions='Medium';
    elseif (strcmp(Selected_SurfaceArea,'Large')==1)&&(strcmp(Selected_NumberOfNodes,'Medium')==1)
        Selected_NumberOfRegions='Medium'; 
    elseif (strcmp(Selected_SurfaceArea,'Large')==1)&&(strcmp(Selected_NumberOfNodes,'Large')==1)
        Selected_NumberOfRegions='Many'; 
    elseif (strcmp(Selected_SurfaceArea,'Large')==1)&&(strcmp(Selected_NumberOfNodes,'VeryLarge')==1)
        Selected_NumberOfRegions='Many';        
    elseif (strcmp(Selected_SurfaceArea,'VeryLarge')==1)&&(strcmp(Selected_NumberOfNodes,'Small')==1)
        Selected_NumberOfRegions='Many';
    elseif (strcmp(Selected_SurfaceArea,'VeryLarge')==1)&&(strcmp(Selected_NumberOfNodes,'Medium')==1)
        Selected_NumberOfRegions='Many'; 
    elseif (strcmp(Selected_SurfaceArea,'VeryLarge')==1)&&(strcmp(Selected_NumberOfNodes,'Large')==1)
        Selected_NumberOfRegions='SoMany'; 
    elseif (strcmp(Selected_SurfaceArea,'VeryLarge')==1)&&(strcmp(Selected_NumberOfNodes,'VeryLarge')==1)
        Selected_NumberOfRegions='SoMany';     
    end;
        
    %          
    if strcmp(Selected_NumberOfRegions,'Feeble')==1
        Region=1;
    elseif strcmp(Selected_NumberOfRegions,'Few')==1
        Region=4;
    elseif strcmp(Selected_NumberOfRegions,'Medium')==1
        Region=9;
    elseif strcmp(Selected_NumberOfRegions,'Many')==1
        Region=16;  
    elseif strcmp(Selected_NumberOfRegions,'SoMany')==1
        Region=25;          
    end;
end

