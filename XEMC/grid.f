c     -read kinematics settings from screen input, use fixed value of dx,dy,dp 
c     -grid on W2 and theta to generate a table that covers the simulation range.
c     -later XEMC  will read the output table and calcuate radiation 
c     correction at each grid point. 
c------------------------------------------------

      program grid
      implicit none

      integer i,stat
      character*5 id,id1
      character*32  output, temp
      real*8 E0, p0, theta_c,dp, yp, xp, w2, q2
      real*8 theta_min,theta_max,tt
      real*8 p_min,p_max,pp,xbj
      logical thend

      real*8 PI, Mp, dw2, dtheta,dE,rad
      parameter(PI=3.14159265, Mp=0.9384)
      parameter(dw2=0.01,dtheta=0.1)
      parameter(rad = PI/180)

      call getarg(1,id)
      write(6,*) "--now generate grid points for run ",id

      thend = .false.

      dp = 10;
      yp = 100;
      xp = 100;

c      write(6,*) "enter E0, p0 in GeV"
c      read(*,*)  E0, p0
c      write(6,*) "enter spec center angle"
c      read(*,*) theta_c


      open(unit=19,file='runplan.inp'
     >                    ,status='old')
      read(19,*)

      dowhile(.not.thend)
        read(19,*) id1, E0,p0,theta_c     
        if(id1==id) then
      
           output = 'table_'//id//'.inp'
           
           
           open(unit=21,iostat=stat,file=output,status='old')
           if(stat==0) close(21,status='delete')
           open(unit=21,file=output,status='new')
           write(21,*) 'E0      p0      theta      Q2      W2   Xbj'
           do i=1,5
              write(21,*) '#'
           enddo
           
           p0 = p0*1000.
           p_min = p0/100. * (100-dp)/1000.
           p_max = p0/100. * (100+dp)/1000.
          

           
          theta_min = acos (cos(theta_c*rad-yp/1000.)*cos(xp/1000.))/rad
          theta_max = acos (cos(theta_c*rad+yp/1000.)*cos(xp/1000.))/rad
          
c          p_min = 3.18; p_max = 3.9;
c          theta_min = 16; theta_max = 27;
          write(6,*) p_max,p_min,theta_max, theta_min           
           tt = theta_c - dtheta*int((theta_c-theta_min)/dtheta+1) !lower limit on theta
           
c           tt = 17
           do while (tt.LE.theta_max)
              
              dE = 2*Mp*E0 - 4*E0*sin(tt*rad/2)**2
              dE = dw2/dE;
              pp = p0-dE*int((p0-p_min)/dE+1)
              
              do while (pp.LE.p_max)
                 q2 = 4*E0*pp * sin(tt*rad/2)**2
                 w2 = Mp**2 + 2*Mp*(E0-pp) - q2
                 xbj = q2/2/0.938/(E0-pp)
c                 if(w2.GT.0) then
                    write(21,211) E0,pp,tt,q2,w2,xbj
c                 endif
                 pp = pp+dE
              enddo
              tt = tt+dtheta
           enddo
           
           
           
           close(21)
           write(6,*) "Done! Generated ", output
 211       format(6(f10.3))
        
        endif
           
        enddo
        
        end program

               




      

