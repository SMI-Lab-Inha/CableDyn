! File: src/CableDyn_FatalReport.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_FatalReport
  !! The standalone driver's last report when the process ends abnormally (src/cabledyn_fatal.c).
  !!
  !! CD_Fatal_Report_Install installs, for the driver process only, the handlers that write
  !! one line to stderr when the run is interrupted (Ctrl+C, a closed console, a termination
  !! signal) or ends in a fatal fault (an access violation, a stack overflow), naming the cause
  !! and the simulated time of the last committed step. The event then continues to the
  !! handler that had it before, so the exit status is unchanged. CD_Fatal_Report_Time records
  !! that time; the time marches call it after every committed step. It only stores a number,
  !! so a host that loads the library without installing the handlers is not affected.
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
  USE CableDyn_Precision, ONLY: wp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CD_Fatal_Report_Install, CD_Fatal_Report_Time

  INTERFACE
    SUBROUTINE c_fatal_report_install() BIND(C, name='cabledyn_fatal_report_install')
    END SUBROUTINE c_fatal_report_install
    SUBROUTINE c_fatal_report_time(simulated_time) BIND(C, name='cabledyn_fatal_report_time')
      IMPORT :: C_DOUBLE
      REAL(C_DOUBLE), VALUE :: simulated_time
    END SUBROUTINE c_fatal_report_time
  END INTERFACE

CONTAINS

  SUBROUTINE CD_Fatal_Report_Install()
    !! Install the abnormal-end report for this process (idempotent). Executables only.
    CALL c_fatal_report_install()
  END SUBROUTINE CD_Fatal_Report_Install

  SUBROUTINE CD_Fatal_Report_Time(simulated_time)
    !! Record the simulated time of the last committed step for the abnormal-end report.
    REAL(wp), INTENT(IN) :: simulated_time
    CALL c_fatal_report_time(REAL(simulated_time, C_DOUBLE))
  END SUBROUTINE CD_Fatal_Report_Time

END MODULE CableDyn_FatalReport
