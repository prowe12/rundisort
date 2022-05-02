

The code disort_driver.f was written for running DISORT from Python using f2py.
Compile as follows:
$ f2py -c -m disort_driver_py disort_driver.f disort_py.f errpack.f linpak.f rdi1mach.f bdref.f --fcompiler=gfortran


I made the following modifications to the standard DISORT code
(these are also contained within the code).
Line numbers are to the left:

1c   Modifications made by P. Rowe:
2c       Parameter MXCLY (max layers) changed from 6 to 120
3c       Commented out the following lines which display header
4c          IF( .NOT.PASS1 .AND. LEN( HEADER ).NE.0 )
5c        &    WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' )
6c        &    ' DISORT: '//HEADER
c
c    Also added msg, errflag to outputs in DISORT call 



Line numbers for changes (may be approximate).

The change to MXCLY was made on line 369.

The commented-out header is on line 493.

I made a variety of changes to add msg and errflag. These are indicated with PMR and are on the following lines:
lines 373-377, 454-457, 511-519, 546-551, 2500-2502, 2545-2550, 2695-2705, 4176-4178, 4213-4218, 4233-4234, 4282-4288, 4922-4932, 5199-5205

