      PROGRAM TEST1
*
*  -- LAPACK test routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Local Scalars ..
      INTEGER            I1, I2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Executable Statements ..
*
*
*     Determine the character set.
*
      I1 = ICHAR( 'A' )
      I2 = ICHAR( 'a' )
      IF( I2-I1.EQ.32 ) THEN
         WRITE( *, * ) ' ASCII character set'
      ELSE
         WRITE( *, * ) ' Non-ASCII character set, IOFF should be ',I2-I1
      END IF
*
*     Test LSAME.
*
      IF( .NOT.LSAME( 'A', 'A' ) )
     $   WRITE( *, 9999 )'A', 'A'
      IF( .NOT.LSAME( 'A', 'a' ) )
     $   WRITE( *, 9999 )'A', 'a'
      IF( .NOT.LSAME( 'a', 'A' ) )
     $   WRITE( *, 9999 )'a', 'A'
      IF( .NOT.LSAME( 'a', 'a' ) )
     $   WRITE( *, 9999 )'a', 'a'
      IF( LSAME( 'A', 'B' ) )
     $   WRITE( *, 9998 )'A', 'B'
      IF( LSAME( 'A', 'b' ) )
     $   WRITE( *, 9998 )'A', 'b'
      IF( LSAME( 'a', 'B' ) )
     $   WRITE( *, 9998 )'a', 'B'
      IF( LSAME( 'a', 'b' ) )
     $   WRITE( *, 9998 )'a', 'b'
      IF( LSAME( 'O', '/' ) )
     $   WRITE( *, 9998 )'O', '/'
      IF( LSAME( '/', 'O' ) )
     $   WRITE( *, 9998 )'/', 'O'
      IF( LSAME( 'o', '/' ) )
     $   WRITE( *, 9998 )'o', '/'
      IF( LSAME( '/', 'o' ) )
     $   WRITE( *, 9998 )'/', 'o'
      WRITE( *, * )' Tests completed'
*
 9999 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is .FALSE.' )
 9998 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is .TRUE.' )
      END
