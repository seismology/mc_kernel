! Martin van Driel, 02/2013
! Martin@vanDriel.de
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This class can be used to print a nice progressbar:
! [==================================================] 100% 
!
! output is only updated if percentage changes (important for performance if
! called very often)

module progressbar_mod
    
    implicit none
    
    private
    type progressbar
        integer     :: last_percentage = -1
        contains
        procedure   :: printbar
    end type
    public :: progressbar, printbar

contains

subroutine printbar(bar, percentage)
    class(progressbar)  :: bar
    integer, intent(in) :: percentage
    integer             :: k

    if (percentage .ne. bar%last_percentage) then
        write(*,'(256a1)', advance='no') (char(8), k =1,100)
        write(*,'("[", 50a1, "] ", i3, "%")', advance='no') &
                ('=', k=1, percentage / 2), &
                (' ', k=1, 50 - percentage / 2), percentage
        bar%last_percentage = percentage
        if (percentage == 100) print *, ''
    endif

end

end module
