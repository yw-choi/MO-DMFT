      module parallel_params

         implicit none
         include 'mpif.h'

         integer, save :: Node
         integer, save :: Nodes
         logical, save :: ParallelOverK

         integer, parameter :: comm = MPI_Comm_world
         integer :: ierr


      end module parallel_params
