#!/bin/csh
#set file = $1
#set wall = $2
make clean;make
#foreach file (00010 00011 00020 00021 00030 00040 00050)
#    foreach i (1 2 3 4 5) 
foreach file(00730 00731 00739 00740 00747 00748 00755 00756 00763)
#foreach file (00730)
    foreach i (6)
    ./mc_hrs_single ar_c ${file} ${i}
    cd worksim
    h2root ar_${file}_${i}.rzdat
    cd ..
    end
       
end
#root ${file}_${wall}.root -l



#
