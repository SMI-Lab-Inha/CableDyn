/* File: tests/fatal_report_faults.c
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
 *
 * Faults that tests/test_fatal_report.f90 provokes to check the abnormal-end report.
 */
#include <signal.h>

void fatal_report_test_null_write(void);
int fatal_report_test_raise_term(void);

static volatile int *volatile fatal_report_null_target = 0;

/* A write through a null pointer: an access violation / SIGSEGV. */
void fatal_report_test_null_write(void)
{
    *fatal_report_null_target = 1;
}

/* SIGTERM to this process; returns 0 where the report does not handle it (Windows). */
int fatal_report_test_raise_term(void)
{
#if defined(_WIN32)
    return 0;
#else
    raise(SIGTERM);
    return 1;
#endif
}
