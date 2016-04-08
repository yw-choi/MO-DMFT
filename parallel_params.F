      module parallel_params

         implicit none
         include 'mpif.h'

         integer, save, public :: Node
         integer, save, public :: Nodes
         logical, save, public :: ParallelOverK

         integer, parameter, public :: comm = MPI_Comm_world
         integer, public :: ierr

         private

      end module parallel_params
