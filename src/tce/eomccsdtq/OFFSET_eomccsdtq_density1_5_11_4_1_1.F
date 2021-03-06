      SUBROUTINE OFFSET_eomccsdtq_density1_5_11_4_1_1(l_a_offset,k_a_off
     &set,size)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i4 ( h1 h2 h7 h10 p5 p6 p9 p13 )_y
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER h1b
      INTEGER h2b
      INTEGER h7b
      INTEGER h10b
      INTEGER p5b
      INTEGER p6b
      INTEGER p9b
      INTEGER p13b
      length = 0
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h7b = h2b,noab
      DO h10b = h7b,noab
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p9b = p6b,noab+nvab
      DO p13b = p9b,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h7b-1)
     &+int_mb(k_spin+h10b-1) .eq. int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b
     &-1)+int_mb(k_spin+p9b-1)+int_mb(k_spin+p13b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_mb(
     &k_sym+h7b-1),ieor(int_mb(k_sym+h10b-1),ieor(int_mb(k_sym+p5b-1),ie
     &or(int_mb(k_sym+p6b-1),ieor(int_mb(k_sym+p9b-1),int_mb(k_sym+p13b-
     &1)))))))) .eq. irrep_y) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1
     &)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h10b-1)+int_mb(k_spin+p5b-1)+
     &int_mb(k_spin+p6b-1)+int_mb(k_spin+p9b-1)+int_mb(k_spin+p13b-1).ne
     &.16)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('eomccsdtq_density1_5_11_4_1_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h7b = h2b,noab
      DO h10b = h7b,noab
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p9b = p6b,noab+nvab
      DO p13b = p9b,noab+nvab
      IF (int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h7b-1)
     &+int_mb(k_spin+h10b-1) .eq. int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b
     &-1)+int_mb(k_spin+p9b-1)+int_mb(k_spin+p13b-1)) THEN
      IF (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),ieor(int_mb(
     &k_sym+h7b-1),ieor(int_mb(k_sym+h10b-1),ieor(int_mb(k_sym+p5b-1),ie
     &or(int_mb(k_sym+p6b-1),ieor(int_mb(k_sym+p9b-1),int_mb(k_sym+p13b-
     &1)))))))) .eq. irrep_y) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1
     &)+int_mb(k_spin+h7b-1)+int_mb(k_spin+h10b-1)+int_mb(k_spin+p5b-1)+
     &int_mb(k_spin+p6b-1)+int_mb(k_spin+p9b-1)+int_mb(k_spin+p13b-1).ne
     &.16)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = p13b - noab - 1 + nvab * (p9b - noab - 1
     & + nvab * (p6b - noab - 1 + nvab * (p5b - noab - 1 + nvab * (h10b 
     &- 1 + noab * (h7b - 1 + noab * (h2b - 1 + noab * (h1b - 1)))))))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_
     &mb(k_range+h7b-1) * int_mb(k_range+h10b-1) * int_mb(k_range+p5b-1)
     & * int_mb(k_range+p6b-1) * int_mb(k_range+p9b-1) * int_mb(k_range+
     &p13b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
