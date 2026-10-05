/* File: src/cabledyn_fatal.c
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
 *
 * A last report from the standalone driver when the process ends abnormally.
 *
 * The driver reports every failure it detects itself on stderr. A process can also end
 * without the program choosing to: an interrupt (Ctrl+C, closing the console window, a
 * termination signal) or a fatal fault (an access violation, a stack overflow). The runtime
 * libraries print nothing for some of these -- a Windows stack overflow in a GNU build ends
 * the process without a word -- and none of them knows how far the simulation got.
 *
 * cabledyn_fatal_report_install() installs handlers that write one line to stderr naming the
 * cause and the simulated time of the last committed step, then let the event continue to
 * whatever handled it before (the Fortran runtime's own report, the operating system's
 * default action), so the exit status is the one the event would have produced anyway.
 * cabledyn_fatal_report_time() records that time; the time march calls it after each
 * committed step. Only the driver installs the handlers: the shared library never changes
 * the signal or exception handling of the host process that loads it.
 *
 * A process that is ended from outside by TerminateProcess (Task Manager, `taskkill /F`) or
 * by SIGKILL (`kill -9`) runs none of its own code, so no program can report that; the
 * driver's closing status line (app/cabledyn.f90) lets a caller recognise such an exit.
 */
/* sigaction(), sigaltstack() and SA_ONSTACK are outside strict C11; request the platform
 * declarations. */
#define _DEFAULT_SOURCE
#define _DARWIN_C_SOURCE
#include <stddef.h>
#include <string.h>

void cabledyn_fatal_report_install(void);
void cabledyn_fatal_report_time(double simulated_time);
double cabledyn_fatal_report_last_time(void);

/* Simulated time of the last committed step; negative until the first one. An aligned
 * double is written in one store on every supported target, and the handlers only read it. */
static volatile double cabledyn_last_time = -1.0;

void cabledyn_fatal_report_time(double simulated_time)
{
    if (simulated_time == simulated_time && simulated_time >= 0.0) {
        cabledyn_last_time = simulated_time;
    }
}

/* The recorded time (negative before the first committed step); for the tests. */
double cabledyn_fatal_report_last_time(void)
{
    return cabledyn_last_time;
}

/* ---- message assembly: no allocation, no stdio, usable from a signal handler ---------- */

typedef struct {
    char text[512];
    size_t len;
} cabledyn_line;

static void line_add(cabledyn_line *line, const char *s)
{
    while (*s != '\0' && line->len + 1 < sizeof line->text) {
        line->text[line->len++] = *s++;
    }
    line->text[line->len] = '\0';
}

static void line_add_unsigned(cabledyn_line *line, unsigned long long value, int min_digits)
{
    char digits[24];
    int n = 0;
    do {
        digits[n++] = (char)('0' + (int)(value % 10ULL));
        value /= 10ULL;
    } while (value != 0ULL && n < (int)sizeof digits);
    while (n < min_digits && n < (int)sizeof digits) {
        digits[n++] = '0';
    }
    while (n > 0 && line->len + 1 < sizeof line->text) {
        line->text[line->len++] = digits[--n];
    }
    line->text[line->len] = '\0';
}

static void line_add_hex(cabledyn_line *line, unsigned long value)
{
    static const char hex[] = "0123456789ABCDEF";
    char digits[2 * sizeof value];
    int n;
    line_add(line, "0x");
    for (n = (int)(2 * sizeof value) - 1; n >= 0; --n) {
        digits[n] = hex[value & 0xFUL];
        value >>= 4;
    }
    for (n = 0; n < (int)(2 * sizeof value) - 8 && digits[n] == '0'; ++n) {
    }
    for (; n < (int)(2 * sizeof value) && line->len + 1 < sizeof line->text; ++n) {
        line->text[line->len++] = digits[n];
    }
    line->text[line->len] = '\0';
}

/* "after the step at simulated time t = 7534.600 s" or the pre-march wording. */
static void line_add_when(cabledyn_line *line)
{
    double t = cabledyn_last_time;
    unsigned long long millis;
    if (!(t >= 0.0) || t > 1.0e15) {
        line_add(line, "before the first time step");
        return;
    }
    millis = (unsigned long long)(t * 1000.0 + 0.5);
    line_add(line, "after the step at simulated time t = ");
    line_add_unsigned(line, millis / 1000ULL, 1);
    line_add(line, ".");
    line_add_unsigned(line, millis % 1000ULL, 3);
    line_add(line, " s");
}

static void line_finish(cabledyn_line *line)
{
    line_add(line, ". The run did not finish; its output files end at the last step written.\n");
}

#if defined(_WIN32)
#include <windows.h>

static LPTOP_LEVEL_EXCEPTION_FILTER cabledyn_previous_filter = NULL;
/* Set by the first report, so one event that reaches two handlers is reported once. */
static volatile LONG cabledyn_reported = 0;

static void write_stderr(const cabledyn_line *line)
{
    HANDLE err = GetStdHandle(STD_ERROR_HANDLE);
    DWORD written = 0;
    if (err != NULL && err != INVALID_HANDLE_VALUE) {
        WriteFile(err, line->text, (DWORD)line->len, &written, NULL);
    }
}

static const char *exception_name(DWORD code)
{
    switch (code) {
    case EXCEPTION_STACK_OVERFLOW:
        return "stack overflow";
    case EXCEPTION_ACCESS_VIOLATION:
        return "access violation";
    case EXCEPTION_IN_PAGE_ERROR:
        return "in-page error";
    case EXCEPTION_ILLEGAL_INSTRUCTION:
        return "illegal instruction";
    case EXCEPTION_PRIV_INSTRUCTION:
        return "privileged instruction";
    case EXCEPTION_INT_DIVIDE_BY_ZERO:
        return "integer division by zero";
    case EXCEPTION_INT_OVERFLOW:
        return "integer overflow";
    case EXCEPTION_FLT_DIVIDE_BY_ZERO:
        return "floating-point division by zero";
    case EXCEPTION_FLT_INVALID_OPERATION:
        return "invalid floating-point operation";
    case EXCEPTION_FLT_OVERFLOW:
        return "floating-point overflow";
    case EXCEPTION_FLT_UNDERFLOW:
        return "floating-point underflow";
    case EXCEPTION_FLT_INEXACT_RESULT:
        return "inexact floating-point result";
    case EXCEPTION_FLT_DENORMAL_OPERAND:
        return "denormal floating-point operand";
    case EXCEPTION_FLT_STACK_CHECK:
        return "floating-point stack check";
    case EXCEPTION_DATATYPE_MISALIGNMENT:
        return "misaligned data access";
    case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        return "array bounds exceeded";
    case EXCEPTION_NONCONTINUABLE_EXCEPTION:
        return "non-continuable exception";
    case 0xC0000409UL: /* STATUS_STACK_BUFFER_OVERRUN: also raised by abort() / fast-fail */
        return "fast-fail abort";
    case 0xC0000374UL: /* STATUS_HEAP_CORRUPTION */
        return "heap corruption";
    default:
        return "unhandled exception";
    }
}

/* Faults no part of the driver recovers from. Floating-point exceptions are not trapped by
 * default; when a build enables traps, the Fortran runtime reports them itself. */
static int is_fatal_fault(DWORD code)
{
    switch (code) {
    case EXCEPTION_STACK_OVERFLOW:
    case EXCEPTION_ACCESS_VIOLATION:
    case EXCEPTION_IN_PAGE_ERROR:
    case EXCEPTION_ILLEGAL_INSTRUCTION:
    case EXCEPTION_PRIV_INSTRUCTION:
    case EXCEPTION_INT_DIVIDE_BY_ZERO:
    case 0xC0000409UL:
    case 0xC0000374UL:
        return 1;
    default:
        return 0;
    }
}

static void report_exception(DWORD code)
{
    if (InterlockedExchange(&cabledyn_reported, 1) == 0) {
        cabledyn_line line;
        line.len = 0;
        line.text[0] = '\0';
        line_add(&line, "\nCableDyn_driver: fatal error: ");
        line_add(&line, exception_name(code));
        line_add(&line, " (exception ");
        line_add_hex(&line, (unsigned long)code);
        line_add(&line, ") ");
        line_add_when(&line);
        line_finish(&line);
        write_stderr(&line);
    }
}

/* First in line for the fatal faults: the GNU runtime handles an access violation in a
 * frame-based handler (its SIGSEGV backtrace) that ends the process before any unhandled-
 * exception filter runs. The search continues afterwards, so the runtime still reports. */
static LONG WINAPI cabledyn_vectored_handler(EXCEPTION_POINTERS *info)
{
    if (info != NULL && info->ExceptionRecord != NULL &&
        is_fatal_fault(info->ExceptionRecord->ExceptionCode)) {
        report_exception(info->ExceptionRecord->ExceptionCode);
    }
    return EXCEPTION_CONTINUE_SEARCH;
}

/* Any other exception nothing handled. */
static LONG WINAPI cabledyn_exception_filter(EXCEPTION_POINTERS *info)
{
    if (info != NULL && info->ExceptionRecord != NULL) {
        report_exception(info->ExceptionRecord->ExceptionCode);
    }
    if (cabledyn_previous_filter != NULL) {
        return cabledyn_previous_filter(info);
    }
    return EXCEPTION_CONTINUE_SEARCH;
}

static BOOL WINAPI cabledyn_console_handler(DWORD event)
{
    const char *cause;
    switch (event) {
    case CTRL_C_EVENT:
        cause = "Ctrl+C";
        break;
    case CTRL_BREAK_EVENT:
        cause = "Ctrl+Break";
        break;
    case CTRL_CLOSE_EVENT:
        cause = "closing its console window";
        break;
    case CTRL_LOGOFF_EVENT:
        cause = "the user logging off";
        break;
    case CTRL_SHUTDOWN_EVENT:
        cause = "the system shutting down";
        break;
    default:
        return FALSE;
    }
    if (InterlockedExchange(&cabledyn_reported, 1) == 0) {
        cabledyn_line line;
        line.len = 0;
        line.text[0] = '\0';
        line_add(&line, "\nCableDyn_driver: stopped by ");
        line_add(&line, cause);
        line_add(&line, " ");
        line_add_when(&line);
        line_finish(&line);
        write_stderr(&line);
    }
    /* Not handled here: the next handler (the Fortran runtime's, or the default one) ends
     * the process with its usual status. */
    return FALSE;
}

void cabledyn_fatal_report_install(void)
{
    static volatile LONG installed = 0;
    ULONG guarantee = 64UL * 1024UL;
    if (InterlockedExchange(&installed, 1) != 0) {
        return;
    }
    /* Keep room on the main thread's stack for the handlers to run after a stack overflow. */
    SetThreadStackGuarantee(&guarantee);
    AddVectoredExceptionHandler(1, cabledyn_vectored_handler);
    cabledyn_previous_filter = SetUnhandledExceptionFilter(cabledyn_exception_filter);
    SetConsoleCtrlHandler(cabledyn_console_handler, TRUE);
}

#else /* POSIX */
#include <signal.h>
#include <unistd.h>

#define CABLEDYN_NSIG 8
static const int cabledyn_signals[CABLEDYN_NSIG] = {SIGSEGV, SIGBUS, SIGFPE, SIGILL,
                                                    SIGABRT, SIGINT, SIGTERM, SIGHUP};
static struct sigaction cabledyn_previous[CABLEDYN_NSIG];
static volatile sig_atomic_t cabledyn_reported = 0;
/* The alternate stack the handlers run on, so a stack overflow can still be reported. */
static char cabledyn_altstack[64 * 1024];

static const char *signal_name(int sig)
{
    switch (sig) {
    case SIGSEGV:
        return "segmentation fault (SIGSEGV; a stack overflow also ends this way)";
    case SIGBUS:
        return "bus error (SIGBUS)";
    case SIGFPE:
        return "floating-point exception (SIGFPE)";
    case SIGILL:
        return "illegal instruction (SIGILL)";
    case SIGABRT:
        return "abort (SIGABRT)";
    case SIGINT:
        return "an interrupt (SIGINT, Ctrl+C)";
    case SIGTERM:
        return "a termination request (SIGTERM)";
    case SIGHUP:
        return "a hang-up (SIGHUP, terminal closed)";
    default:
        return "a signal";
    }
}

static void cabledyn_signal_handler(int sig)
{
    int k;
    if (!cabledyn_reported) {
        cabledyn_line line;
        int fault = sig == SIGSEGV || sig == SIGBUS || sig == SIGFPE || sig == SIGILL ||
                    sig == SIGABRT;
        cabledyn_reported = 1;
        line.len = 0;
        line.text[0] = '\0';
        line_add(&line, fault ? "\nCableDyn_driver: fatal error: " : "\nCableDyn_driver: stopped by ");
        line_add(&line, signal_name(sig));
        line_add(&line, " ");
        line_add_when(&line);
        line_finish(&line);
        {
            ssize_t ignored = write(STDERR_FILENO, line.text, line.len);
            (void)ignored;
        }
    }
    /* Hand the signal to whatever handled it before (the Fortran runtime's backtrace or the
     * default action), so the exit status is unchanged. SA_NODEFER delivers it at once. */
    for (k = 0; k < CABLEDYN_NSIG; ++k) {
        if (cabledyn_signals[k] == sig) {
            sigaction(sig, &cabledyn_previous[k], NULL);
            break;
        }
    }
    raise(sig);
}

void cabledyn_fatal_report_install(void)
{
    static int installed = 0;
    stack_t alt;
    struct sigaction action;
    int k;
    if (installed) {
        return;
    }
    installed = 1;
    memset(&alt, 0, sizeof alt);
    alt.ss_sp = cabledyn_altstack;
    alt.ss_size = sizeof cabledyn_altstack;
    alt.ss_flags = 0;
    sigaltstack(&alt, NULL);
    memset(&action, 0, sizeof action);
    action.sa_handler = cabledyn_signal_handler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_ONSTACK | SA_NODEFER;
    for (k = 0; k < CABLEDYN_NSIG; ++k) {
        if (sigaction(cabledyn_signals[k], NULL, &cabledyn_previous[k]) != 0) {
            continue;
        }
        /* A signal the parent set to be ignored (nohup, background job) stays ignored. */
        if (cabledyn_previous[k].sa_handler == SIG_IGN) {
            continue;
        }
        sigaction(cabledyn_signals[k], &action, NULL);
    }
}
#endif
