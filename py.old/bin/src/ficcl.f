        SUBROUTINE CLOSING
        INCLUDE 'ficc.f'
!
!
        CLOSE (2)
        CLOSE (11)
        CLOSE (12)
        CLOSE (13)
        CLOSE (15)
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
        CLOSE (19)
        CLOSE (41)
!       CLOSE (70)
!
        CLOSE (3)
        CLOSE (10)
        CLOSE (20)
        CLOSE (21)
        CLOSE (29)
        CLOSE (30)
        CLOSE (31)
        CLOSE (32)
        CLOSE (42)
        CLOSE (43)
        CLOSE (46)
!       SHChange: Close additional files	
        IF (WHIST==1) CLOSE (50)
!        CLOSE (51)
!        CLOSE (52)
!        CLOSE (53)
!        CLOSE (54)
!        CLOSE (55)
!        CLOSE (56)
!        CLOSE (57)
!        CLOSE (58)
!        CLOSE (61)
!        CLOSE (62)
!        CLOSE (63)
!        CLOSE (64)
!        CLOSE (65)
!        CLOSE (66)
!        CLOSE (67)
!        CLOSE (68)
!        CLOSE (69)

!
        RETURN
        END     
                
                
