program testuniform

      real, dimension(3,4)   :: vertices

      real, dimension(3, 100) :: points



      do
         do ivert = 1,4
            print *, 'Enter Coordinates of point ', ivert
            read(*,*) vertices(:,ivert)
         end do


         call uniform_in_tetrahedron ( vertices, 100, 123, points )

         do ipoint = 1,100
            print *, points(:,ipoint)
         end do
      end do

      end program
