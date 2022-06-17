      subroutine test(n)
      integer n

c     Lines added for f2py:
cf2py intent(in)  :: n

        do 100 i=0, n
            print *, "Fortran says hello"
100     continue
        end