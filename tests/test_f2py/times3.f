      subroutine times3(n, m)
      integer n, m

c     Lines added for f2py:
cf2py intent(in)  :: n
cf2py intent(out) :: m

        m = n * 3

        end