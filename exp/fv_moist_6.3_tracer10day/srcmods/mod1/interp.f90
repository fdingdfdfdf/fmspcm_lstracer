subroutine findind(nx, xa, x, ind)

implicit none

integer, intent(in) :: nx
real(8), intent(in), dimension(nx) :: xa
real(8), intent(in) :: x
integer, intent(out) :: ind

integer :: i
real :: diff, difft

ind = 1
diff = abs(xa(ind) - x)

do i = 1, nx
   difft = abs(xa(i) - x)
   if (difft < diff) then
      diff = difft
      ind = i
   end if
end do

end subroutine findind

subroutine interp(nx, xa, ya, mx, xx, yy)
! xa is increasing
implicit none

integer, intent(in) :: nx, mx
real(8), intent(in), dimension(nx) :: xa, ya
real(8), intent(in), dimension(mx) :: xx
real(8), intent(out), dimension(mx) :: yy

integer :: nn = 4
integer :: i1, i2, j, n, ns, m, i, ind
real(8), allocatable, dimension(:) :: b, c, d
real(8) :: diff, difft, x, y, ho, hp, w, dy

do j = 1, mx

   x = xx(j)

   call findind(nx, xa, x, ind)

   i1 = max(ind - nn, 1)
   i2 = min(ind + nn, nx)

   n = i2 - i1 + 1

   allocate(b(n))
   allocate(c(n))
   allocate(d(n))

   b = xa(i1:i2)
   c = ya(i1:i2)
   d = ya(i1:i2)

   ns = 1
   diff = abs(b(1) - x)

   do i = 2, n
      difft = abs(b(i) - x)
      if (difft < diff) then
         diff = difft
         ns = i
      end if
   end do

   y = c(ns)
   ns = ns - 1
   do m = 1, n-1
      do i = 1, (n-m)
         ho = b(i) - x
         hp = b(i+m) - x
         w = c(i+1) - d(i)
         c(i) = ho * w/(ho - hp)
         d(i) = hp * w/(ho - hp)
      end do
      
      if (2*ns .lt. (n-m)) then
         dy = c(ns+1)
      else
         dy = d(ns)
         ns = ns - 1
      end if
      y = y + dy
   end do
   deallocate(b, c, d)

   yy(j) = y

end do

end subroutine interp
