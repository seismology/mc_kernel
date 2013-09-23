module ab
  implicit none
  public
contains

integer function a_plus_b(a,b)
  integer, intent(in)   :: a, b

  a_plus_b = a + b
end function

end module
