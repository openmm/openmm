c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine calendar  --  find the current date and time  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "calendar" returns the current time as a set of integer values
c     representing the year, month, day, hour, minute and second
c
c     note only one of the various implementations below should
c     be activated by removing comment characters
c
c
      subroutine calendar (year,month,day,hour,minute,second)
      implicit none
      integer year,month
      integer day,hour
      integer minute,second
c
c
c     use the standard "date_and_time" intrinsic function
c
      integer values(8)
      character*5 zone
      character*8 date
      character*10 time
      call date_and_time (date,time,zone,values)
      year = values(1)
      month = values(2)
      day = values(3)
      hour = values(5)
      minute = values(6)
      second = values(7)
c
c     use the obsolete "itime" and "idate" intrinsic functions
c
c     integer hms(3)
c     call itime (hms)
c     hour = hms(1)
c     minute = hms(2)
c     second = hms(3)
c     call idate (month,day,year)
      return
      end
