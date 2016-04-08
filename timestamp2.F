      SUBROUTINE timestamp2(str)


#ifdef MPI
#ifdef MPI0
      implicit none
      include 'mpif.h'
#else
      use mpi
      implicit none
#endif
      integer
     .  MPIerror
#endif 


      CHARACTER(len=*) str

      character(len=1), parameter :: dash = "-"
      character(len=1), parameter :: colon = ":"
      character(len=3), parameter :: prefix = ">> "
      character(len=3)  :: month_str(12)
      data month_str 
     $    /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &     'OCT','NOV','DEC'/

      integer sec, min, hour, day, month, year
      integer values(8),Node,Nodes

#ifdef MPI       
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else  
      Node = 0
      Nodes = 1
#endif 

#ifdef MPI
      call MPI_barrier(MPI_Comm_World,MPIerror)
#endif
      if(Node.eq.0) then
         call date_and_time(values=values)
         year = values(1)
         month = values(2)
         day = values(3)
         hour = values(5)
         min = values(6)
         sec = values(7)
         write(6,1000) prefix, trim(str), colon,
     $                 day, dash, month_str(month), dash, year,
     $                 hour, colon, min, colon, sec

 1000    format(2a,a1,2x,i2,a1,a3,a1,i4,2x,i2,a1,i2.2,a1,i2.2)
         call flush_(6)           
      endif
#ifdef MPI
      call MPI_barrier(MPI_Comm_World,MPIerror)
#endif
      END

      SUBROUTINE printi(str,i)


#ifdef MPI
#ifdef MPI0
      implicit none
      include 'mpif.h'
#else
      use mpi
      implicit none
#endif
      integer
     .  MPIerror
#endif 


      CHARACTER(len=*) str

      integer i,Node,Nodes

#ifdef MPI       
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else  
      Node = 0
      Nodes = 1
#endif 

      write(6,'(a,i3,a,a,i10)')' Node = ',Node,'  ',str,i
      call flush_(6)           
      END

