/* File: tests/fatal_report_faults.c
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
 *
 * Faults that tests/test_fatal_report.f90 provokes to check the abnormal-end report.
 */
#define _DEFAULT_SOURCE
#define _DARWIN_C_SOURCE
#include <signal.h>

void fatal_report_test_null_write(void);
int fatal_report_test_raise_term(void);
int fatal_report_test_fault_context(void);
int fatal_report_test_dispositions(void);

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

#if defined(_WIN32)
/* The fault-context check is POSIX-only: returns -1 (skipped). */
int fatal_report_test_fault_context(void)
{
    return -1;
}

/* The disposition check is POSIX-only: returns -1 (skipped). */
int fatal_report_test_dispositions(void)
{
    return -1;
}
#else
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

void cabledyn_fatal_report_install(void);
void cabledyn_fatal_report_time(double simulated_time);

#define FAULT_ADDRESS ((uintptr_t)16)

typedef struct {
    int sig;
    int code;
    uintptr_t addr;
} fault_record;

static int fault_pipe = -1;
/* Read through a volatile, so the compiler cannot treat the store as a known-invalid one. */
static volatile uintptr_t fault_address = FAULT_ADDRESS;

/* A handler installed before the report, as a runtime's would be: it records the siginfo it
 * receives, restores the default action and returns, so the fault re-executes and ends the
 * process. */
static void prior_handler(int sig, siginfo_t *info, void *context)
{
    fault_record record;
    ssize_t ignored;
    (void)context;
    record.sig = sig;
    record.code = info != NULL ? info->si_code : -99;
    record.addr = info != NULL ? (uintptr_t)info->si_addr : 0;
    ignored = write(fault_pipe, &record, sizeof record);
    (void)ignored;
    signal(sig, SIG_DFL);
}

static void read_all(int fd, char *buffer, size_t size)
{
    size_t used = 0;
    ssize_t got;
    while (used + 1 < size && (got = read(fd, buffer + used, size - 1 - used)) > 0) {
        used += (size_t)got;
    }
    buffer[used] = '\0';
}

static int count(const char *text, const char *what)
{
    int n = 0;
    const char *at = text;
    while ((at = strstr(at, what)) != NULL) {
        ++n;
        at += strlen(what);
    }
    return n;
}

/* A child process with a three-argument SIGSEGV/SIGBUS handler installs the report and writes
 * to address 16. The previous handler must receive the hardware fault's own siginfo (a
 * positive si_code and si_addr 16, not a re-sent signal), the child must end by that signal,
 * and the report must be written once. Returns 0 when all hold, 1 otherwise. */
int fatal_report_test_fault_context(void)
{
    int records[2], errors[2], status = 0, ok;
    pid_t child;
    fault_record record;
    ssize_t got;
    char text[8192];
    struct sigaction action;

    if (pipe(records) != 0 || pipe(errors) != 0) {
        fprintf(stderr, "FAIL: fault context: pipe() failed\n");
        return 1;
    }
    fflush(NULL);
    child = fork();
    if (child < 0) {
        fprintf(stderr, "FAIL: fault context: fork() failed\n");
        return 1;
    }
    if (child == 0) {
        close(records[0]);
        close(errors[0]);
        dup2(errors[1], STDERR_FILENO);
        fault_pipe = records[1];
        memset(&action, 0, sizeof action);
        action.sa_sigaction = prior_handler;
        sigemptyset(&action.sa_mask);
        action.sa_flags = SA_SIGINFO;
        sigaction(SIGSEGV, &action, NULL);
        sigaction(SIGBUS, &action, NULL);
        cabledyn_fatal_report_install();
        cabledyn_fatal_report_time(7534.6);
        *(volatile int *)fault_address = 1;
        _exit(3);
    }
    close(records[1]);
    close(errors[1]);
    read_all(errors[0], text, sizeof text);
    got = read(records[0], &record, sizeof record);
    waitpid(child, &status, 0);
    fprintf(stderr, "%s", text);
    if (got != (ssize_t)sizeof record) {
        fprintf(stderr, "FAIL: fault context: the previous handler was not called\n");
        return 1;
    }
    fprintf(stderr,
            "fault context: previous handler got signal %d, si_code %d, si_addr %#lx; child "
            "ended %s %d\n",
            record.sig, record.code, (unsigned long)record.addr,
            WIFSIGNALED(status) ? "by signal" : "with status",
            WIFSIGNALED(status) ? WTERMSIG(status) : WEXITSTATUS(status));
    ok = (record.sig == SIGSEGV || record.sig == SIGBUS) && record.code > 0 &&
         record.addr == FAULT_ADDRESS && WIFSIGNALED(status) && WTERMSIG(status) == record.sig &&
         count(text, "CableDyn_driver: fatal error:") == 1 && count(text, "t = 7534.600 s") == 1;
    if (!ok) {
        fprintf(stderr, "FAIL: fault context: the fault did not keep its own context\n");
        return 1;
    }
    return 0;
}

/* SIGTERM handler of the "handler" case: a real one-argument function. */
static void prior_term_handler(int sig)
{
    static const char text[] = "prior SIGTERM handler ran\n";
    ssize_t ignored = write(STDERR_FILENO, text, sizeof text - 1);
    (void)ignored;
    (void)sig;
    _exit(7);
}

/* Run one child that sets SIGTERM to `setup`, installs the report and sends itself SIGTERM.
 * setup 0: SIG_IGN with SA_SIGINFO; 1: SIG_IGN; 2: SIG_DFL; 3: a one-argument handler. */
static int disposition_case(int setup, const char *name)
{
    int errors[2], status = 0, ok, reports, ran;
    pid_t child;
    char text[4096];
    struct sigaction action;

    if (pipe(errors) != 0) {
        fprintf(stderr, "FAIL: dispositions: pipe() failed\n");
        return 1;
    }
    fflush(NULL);
    child = fork();
    if (child < 0) {
        fprintf(stderr, "FAIL: dispositions: fork() failed\n");
        return 1;
    }
    if (child == 0) {
        close(errors[0]);
        dup2(errors[1], STDERR_FILENO);
        memset(&action, 0, sizeof action);
        sigemptyset(&action.sa_mask);
        if (setup == 0) {
            action.sa_handler = SIG_IGN; /* the sentinel, stored with SA_SIGINFO set */
            action.sa_flags = SA_SIGINFO;
        } else if (setup == 1) {
            action.sa_handler = SIG_IGN;
        } else if (setup == 2) {
            action.sa_handler = SIG_DFL;
        } else {
            action.sa_handler = prior_term_handler;
        }
        sigaction(SIGTERM, &action, NULL);
        cabledyn_fatal_report_install();
        cabledyn_fatal_report_time(7534.6);
        kill(getpid(), SIGTERM);
        /* Still running: an ignored SIGTERM must leave the process alone. */
        _exit(0);
    }
    close(errors[1]);
    read_all(errors[0], text, sizeof text);
    waitpid(child, &status, 0);
    reports = count(text, "CableDyn_driver: stopped by a termination request");
    ran = count(text, "prior SIGTERM handler ran");
    if (setup <= 1) {
        ok = WIFEXITED(status) && WEXITSTATUS(status) == 0 && reports == 0;
    } else if (setup == 2) {
        ok = WIFSIGNALED(status) && WTERMSIG(status) == SIGTERM && reports == 1;
    } else {
        ok = WIFEXITED(status) && WEXITSTATUS(status) == 7 && reports == 1 && ran == 1;
    }
    fprintf(stderr, "dispositions: %s: child ended %s %d, %d report(s)%s\n", name,
            WIFSIGNALED(status) ? "by signal" : "with status",
            WIFSIGNALED(status) ? WTERMSIG(status) : WEXITSTATUS(status), reports,
            ok ? "" : " -- unexpected");
    if (!ok) {
        fprintf(stderr, "%s", text);
    }
    return ok ? 0 : 1;
}

/* An ignored SIGTERM (with or without SA_SIGINFO) stays ignored and is never reported or
 * called; a default one is reported and ends the process by SIGTERM; a previous handler is
 * reported and then runs. Returns 0 when all hold, 1 otherwise. */
int fatal_report_test_dispositions(void)
{
    int failures = 0;
    failures += disposition_case(0, "SIG_IGN with SA_SIGINFO");
    failures += disposition_case(1, "SIG_IGN");
    failures += disposition_case(2, "SIG_DFL");
    failures += disposition_case(3, "a one-argument handler");
    if (failures != 0) {
        fprintf(stderr, "FAIL: dispositions: %d case(s) wrong\n", failures);
        return 1;
    }
    return 0;
}
#endif
