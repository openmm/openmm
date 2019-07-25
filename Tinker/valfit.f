c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module valfit  --  valence term parameter fitting values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     fit_bond    logical flag to fit bond stretch parameters
c     fit_angle   logical flag to fit angle bend parameters
c     fit_strbnd  logical flag to fit stretch-bend parameters
c     fit_urey    logical flag to fit Urey-Bradley parameters
c     fit_opbend  logical flag to fit out-of-plane bend parameters
c     fit_tors    logical flag to fit torsional parameters
c     fit_struct  logical flag to structure-fit valence parameters
c     fit_force   logical flag to force-fit valence parameters
c
c
      module valfit
      implicit none
      logical fit_bond,fit_angle
      logical fit_strbnd,fit_urey
      logical fit_opbend,fit_tors
      logical fit_struct,fit_force
      save
      end
