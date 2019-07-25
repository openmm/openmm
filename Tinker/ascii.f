c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module ascii  --  selected ASCII character code values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     null         decimal value of ASCII code for null (0)
c     tab          decimal value of ASCII code for tab (9)
c     linefeed     decimal value of ASCII code for linefeed (10)
c     formfeed     decimal value of ASCII code for formfeed (12)
c     carriage     decimal value of ASCII code for carriage return (13)
c     escape       decimal value of ASCII code for escape (27)
c     space        decimal value of ASCII code for blank space (32)
c     exclamation  decimal value of ASCII code for exclamation (33)
c     quote        decimal value of ASCII code for double quote (34)
c     pound        decimal value of ASCII code for pound sign (35)
c     dollar       decimal value of ASCII code for dollar sign (36)
c     percent      decimal value of ASCII code for percent sign (37)
c     ampersand    decimal value of ASCII code for ampersand (38)
c     apostrophe   decimal value of ASCII code for single quote (39)
c     asterisk     decimal value of ASCII code for asterisk (42)
c     plus         decimal value of ASCII code for plus sign (43)
c     comma        decimal value of ASCII code for comma (44)
c     minus        decimal value of ASCII code for minus sign (45)
c     period       decimal value of ASCII code for period (46)
c     frontslash   decimal value of ASCII codd for frontslash (47)
c     colon        decimal value of ASCII code for colon (58)
c     semicolon    decimal value of ASCII code for semicolon (59)
c     equal        decimal value of ASCII code for equal sign (61)
c     question     decimal value of ASCII code for question mark (63)
c     atsign       decimal value of ASCII code for at sign (64)
c     backslash    decimal value of ASCII code for backslash (92)
c     caret        decimal value of ASCII code for caret (94)
c     underbar     decimal value of ASCII code for underbar (95)
c     vertical     decimal value of ASCII code for vertical bar (124)
c     tilde        decimal value of ASCII code for tilde (126)
c
c
      module ascii
      implicit none
      integer null,tab
      integer linefeed,formfeed
      integer carriage,escape
      integer space,exclamation
      integer quote,pound
      integer dollar,percent
      integer ampersand
      integer apostrophe
      integer asterisk,plus
      integer comma,minus
      integer period,frontslash
      integer colon,semicolon
      integer equal,question
      integer atsign,backslash
      integer caret,underbar
      integer vertical,tilde
      parameter (null=0)
      parameter (tab=9)
      parameter (linefeed=10)
      parameter (formfeed=12)
      parameter (carriage=13)
      parameter (escape=27)
      parameter (space=32)
      parameter (exclamation=33)
      parameter (quote=34)
      parameter (pound=35)
      parameter (dollar=36)
      parameter (percent=37)
      parameter (ampersand=38)
      parameter (apostrophe=39)
      parameter (asterisk=42)
      parameter (plus=43)
      parameter (comma=44)
      parameter (minus=45)
      parameter (period=46)
      parameter (frontslash=47)
      parameter (colon=58)
      parameter (semicolon=59)
      parameter (equal=61)
      parameter (question=63)
      parameter (atsign=64)
      parameter (backslash=92)
      parameter (caret=94)
      parameter (underbar=95)
      parameter (vertical=124)
      parameter (tilde=126)
      save
      end
