      program test
      integer n

c     Lines added for f2py:
cf2py intent(in)  :: n

        n = 3
        do 100 i=0, n
            print *, "Fortran says hello"
100     continue
        end