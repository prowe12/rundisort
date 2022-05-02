

I made the following modifications to the standard DISORT code
(these are also contained within the code).
Line numbers are to the left:

1c   Modifications made by P. Rowe:
2c       Parameter MXCLY (max layers) changed from 6 to 120
3c       Commented out the following lines which display header
4c          IF( .NOT.PASS1 .AND. LEN( HEADER ).NE.0 )
5c        &    WRITE( *,'(//,1X,100(''*''),/,A,/,1X,100(''*''))' )
6c        &    ' DISORT: '//HEADER

4235c                ** Below, I change the if statement to an if then statement
4236c                ** PMR 2020/01/04

