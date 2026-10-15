/* File: tests/test_c_api_host_coupling.c
 * SPDX-License-Identifier: Apache-2.0
 * Host-coupling gate of C ABI 1 minor extension 2, through the public header
 * and the shared library.
 *
 * Usage: test_c_api_host_coupling <box_points.dat> <mixed_chain_cable.dat>
 *
 * box_points.dat is the foamMooring EsfLOWC box (four slack EI=0 lines, point
 * route) whose fairleads sit at z = -0.062 m in the deck while the CFD box
 * starts with them at z = -0.077 m. mixed_chain_cable.dat is an EI=0 chain
 * plus a finite-EI cable (aggregate route, staged steps). The test checks:
 * - CableDyn_InitDeckAt solves the statics at the host positions: on both
 *   routes the result equals CableDyn_InitDeck of the deck moved there bit for
 *   bit; for the box, line 4's fairlead load is the standalone value and held
 *   fairleads stay at rest, while InitDeck at the deck positions followed by
 *   steps at the host positions rings;
 * - CableDyn_SaveState/RestoreState: a restore followed by steps reproduces a
 *   fresh model taking the same steps bit for bit, on both routes, also after
 *   several steps between save and restore and after repeated restores;
 * - CableDyn_GetCoupledLayout reports anchors and fairleads on the point route
 *   and only the Coupled points on the aggregate route;
 * - the fail-closed rejections of the new calls.
 */
#include "CableDyn_CAPI.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition, message)                        \
    do {                                                 \
        if (!(condition)) {                              \
            fprintf(stderr, "FAIL: %s\n", (message));    \
            return 1;                                    \
        }                                                \
    } while (0)

#define BOX_NDOF 24  /* 8 line-bound points: 4 Fixed anchors, 4 Coupled fairleads */
#define BOX_NHOST 12 /* the 4 Coupled fairleads */
#define BOX_DT 1.0e-3
#define BOX_NSTEP 1000 /* 1 s */
#define DECK_Z (-0.062)
#define HOST_Z (-0.077)

#define MIX_NDOF 6 /* the chain and the cable fairleads; the anchors are not coupled */
#define MIX_DT 0.01

static const char *g_box;
static const char *g_mixed;

static void report(void *h, const char *where, int err)
{
    char msg[1025];
    CableDyn_GetLastError(h, msg, (int)sizeof msg);
    fprintf(stderr, "FAIL: %s (status %d): %s\n", where, err, msg);
}

static int init_deck(void **h, const char *deck)
{
    int err = -1;
    CableDyn_Create(h, &err);
    CHECK(err == CD_C_OK, "create");
    CableDyn_InitDeck(*h, deck, (int)strlen(deck), &err);
    if (err != CD_C_OK) {
        report(*h, "InitDeck", err);
        return 1;
    }
    return 0;
}

static int init_deck_at(void **h, const char *deck, const double *q0, int n_q0)
{
    int err = -1;
    CableDyn_Create(h, &err);
    CHECK(err == CD_C_OK, "create");
    CableDyn_InitDeckAt(*h, deck, (int)strlen(deck), q0, n_q0, &err);
    if (err != CD_C_OK) {
        report(*h, "InitDeckAt", err);
        return 1;
    }
    return 0;
}

static int step(void *h, double dt, const double *q, const double *v, const double *a, int n)
{
    bool converged = false, stalled = true;
    int n_iter = 0, err = -1;
    CableDyn_Step(h, dt, q, v, a, n, &converged, &stalled, &n_iter, &err);
    if (err != CD_C_OK || !converged) {
        report(h, "Step", err);
        return 1;
    }
    return 0;
}

static int loads(void *h, double *f, int n)
{
    int err = -1;
    CableDyn_CalcOutput(h, f, n, &err);
    if (err != CD_C_OK) {
        report(h, "CalcOutput", err);
        return 1;
    }
    return 0;
}

/* The box fairlead positions (deck x, y) at height z, as the 4 Coupled blocks. */
static void box_fairleads(double z, double *q_host)
{
    static const double xy[4][2] = {{-0.1, 0.1}, {-0.1, -0.1}, {0.1, 0.1}, {0.1, -0.1}};
    for (int p = 0; p < 4; ++p) {
        q_host[3 * p + 0] = xy[p][0];
        q_host[3 * p + 1] = xy[p][1];
        q_host[3 * p + 2] = z;
    }
}

/* Copy a deck with the Coupled points moved (the standalone reference of InitDeckAt): on every
 * line naming a Coupled point, the text old_text is replaced by new_text of the same length;
 * n_moved such lines must be found. */
static int write_moved_deck(const char *in_path, const char *out_path, const char *old_text,
                            const char *new_text, int n_moved)
{
    FILE *in = fopen(in_path, "r"), *out = NULL;
    char line[512];
    int moved = 0;
    CHECK(in != NULL && strlen(old_text) == strlen(new_text), "open the deck to move");
    out = fopen(out_path, "w");
    if (out == NULL) {
        fclose(in);
        CHECK(0, "write the moved deck");
    }
    while (fgets(line, (int)sizeof line, in) != NULL) {
        char *at = strstr(line, old_text);
        if (at != NULL && strstr(line, "Coupled") != NULL) {
            memcpy(at, new_text, strlen(new_text));
            ++moved;
        }
        fputs(line, out);
    }
    fclose(in);
    fclose(out);
    CHECK(moved == n_moved, "the deck names its Coupled points as expected");
    return 0;
}

static int check_version(void)
{
    CHECK(CABLEDYN_CAPI_ABI_MINOR == 2 && CableDyn_GetAbiMinor() == CABLEDYN_CAPI_ABI_MINOR,
          "ABI minor level 2");
    return 0;
}

static int check_point_layout(void)
{
    void *h = NULL;
    int err = -1, ids[8], types[8], n;
    if (init_deck(&h, g_box)) return 1;
    n = CableDyn_NCoupledDOF(h, &err);
    CHECK(err == CD_C_OK && n == BOX_NDOF, "the box coupled vector has 8 blocks");
    CableDyn_GetCoupledLayout(h, ids, types, n / 3, &err);
    CHECK(err == CD_C_OK, "layout query");
    for (int k = 0; k < 8; ++k) {
        CHECK(ids[k] == k + 1, "point-route blocks follow the POINTS table");
        CHECK(types[k] == (k < 4 ? CD_C_POINT_FIXED : CD_C_POINT_COUPLED),
              "point-route blocks: 4 Fixed anchors, then 4 Coupled fairleads");
    }
    CableDyn_GetCoupledLayout(h, NULL, types, n / 3, &err);
    CHECK(err == CD_C_OK && types[7] == CD_C_POINT_COUPLED, "a NULL output is skipped");
    CableDyn_GetCoupledLayout(h, ids, types, n / 3 - 1, &err);
    CHECK(err == CD_C_BAD_INPUT, "a wrong block count is rejected");
    CableDyn_Close(&h, &err);
    return 0;
}

static int check_raw_layout(void)
{
    void *h = NULL;
    int err = -1, ids[2], types[2], bogus = 0;
    const double q0[6] = {0.0, 0.0, 0.0, 10.0, 0.0, 0.0}, v0[6] = {0.0};
    const int elem_conn[2] = {1, 2}, fixed[6] = {1, 2, 3, 4, 5, 6};
    const double l0[1] = {10.0}, ea[1] = {1.0e6}, rho_a[1] = {10.0};
    char msg[256];

    CableDyn_Create(&h, &err);
    CableDyn_GetCoupledLayout(h, ids, types, 2, &err);
    CHECK(err == CD_C_NOT_INITIALIZED, "layout of an uninitialized handle");
    CableDyn_SaveState(h, &err);
    CHECK(err == CD_C_NOT_INITIALIZED, "save on an uninitialized handle");
    CableDyn_RestoreState(h, &err);
    CHECK(err == CD_C_NOT_INITIALIZED, "restore on an uninitialized handle");
    CableDyn_InitLine(h, 2, 1, q0, v0, elem_conn, l0, ea, rho_a, fixed, 6, 1, 0.8, &err);
    CHECK(err == CD_C_OK, "raw line init");
    CableDyn_GetCoupledLayout(h, ids, types, 2, &err);
    CableDyn_GetLastError(h, msg, (int)sizeof msg);
    CHECK(err == CD_C_BAD_INPUT && strstr(msg, "raw arrays") != NULL, "a raw-array handle has no layout");
    CableDyn_GetCoupledLayout((void *)&bogus, ids, types, 2, &err);
    CHECK(err == CD_C_BAD_HANDLE, "foreign handle");
    CableDyn_Close(&h, &err);
    return 0;
}

static int check_init_at(void)
{
    void *at = NULL, *full = NULL, *moved = NULL, *plain = NULL;
    double q_host[BOX_NHOST], q_full[BOX_NDOF], q[BOX_NDOF], v[BOX_NDOF] = {0.0}, a[BOX_NDOF] = {0.0};
    double f_at[BOX_NDOF], f_full[BOX_NDOF], f_moved[BOX_NDOF], f[BOX_NDOF], f_plain0[BOX_NDOF];
    double drift = 0.0, ring = 0.0, ring_end = 0.0;
    int err = -1, n, ids[8], types[8];
    const char *moved_deck = "capi_host_box_moved.dat";
    /* standalone CableDyn with the deck moved to z = -0.077 (MoorDyn: 0.3034, -0.0763, -0.4136) */
    const double f4_ref[3] = {0.3026, -0.0761, -0.4114};
    char msg[1025];

    box_fairleads(HOST_Z, q_host);
    if (init_deck_at(&at, g_box, q_host, BOX_NHOST)) return 1;
    n = CableDyn_NCoupledDOF(at, &err);
    CHECK(err == CD_C_OK && n == BOX_NDOF, "InitDeckAt keeps the route's coupled vector");
    CableDyn_GetCoupledLayout(at, ids, types, 8, &err);
    CHECK(err == CD_C_OK && ids[7] == 8 && types[7] == CD_C_POINT_COUPLED, "layout after InitDeckAt");
    CableDyn_GetCoupledMotion(at, q, v, a, n, &err);
    CHECK(err == CD_C_OK, "coupled motion after InitDeckAt");
    for (int p = 0; p < 4; ++p) {
        CHECK(q[3 * p + 2] == -0.5, "the anchors stay at the deck coordinates");
        CHECK(memcmp(&q[12 + 3 * p], &q_host[3 * p], 3 * sizeof(double)) == 0,
              "the fairleads start exactly at q0");
    }
    for (int i = 0; i < n; ++i) CHECK(v[i] == 0.0 && a[i] == 0.0, "the statics start at rest");
    if (loads(at, f_at, n)) return 1;
    printf("InitDeckAt z=%.3f: line 4 fairlead load (%.4f, %.4f, %.4f) N\n", HOST_Z, f_at[21], f_at[22],
           f_at[23]);
    for (int c = 0; c < 3; ++c)
        CHECK(fabs(f_at[21 + c] - f4_ref[c]) < 1.0e-4, "line 4 fairlead load matches the standalone statics");

    /* The full-vector form: the Fixed blocks are not read. */
    memcpy(q_full, q, sizeof q_full);
    for (int i = 0; i < 12; ++i) q_full[i] = 1.0e3;
    if (init_deck_at(&full, g_box, q_full, BOX_NDOF)) return 1;
    if (loads(full, f_full, n)) return 1;
    CHECK(memcmp(f_full, f_at, sizeof f_at) == 0, "full and Coupled-only q0 give the same statics");
    CableDyn_GetCoupledMotion(full, q_full, v, a, n, &err);
    CHECK(err == CD_C_OK && memcmp(q_full, q, sizeof q) == 0, "the full form ignores the Fixed blocks");

    /* Bit for bit the statics of the deck moved to the host positions. */
    if (write_moved_deck(g_box, moved_deck, "-0.062", "-0.077", 4)) return 1;
    if (init_deck(&moved, moved_deck)) return 1;
    if (loads(moved, f_moved, n)) return 1;
    CHECK(memcmp(f_moved, f_at, sizeof f_at) == 0, "InitDeckAt equals InitDeck of the moved deck");
    remove(moved_deck);

    /* Held at the host positions, the InitDeckAt model stays at equilibrium for 1 s. */
    memset(v, 0, sizeof v);
    memset(a, 0, sizeof a);
    for (int s = 1; s <= BOX_NSTEP; ++s) {
        if (step(at, BOX_DT, q, v, a, n) || loads(at, f, n)) return 1;
        for (int i = 12; i < n; ++i) drift = fmax(drift, fabs(f[i] - f_at[i]));
    }
    printf("InitDeckAt, fairleads held 1 s: max fairlead load change %.3e N\n", drift);
    CHECK(drift < 1.0e-6, "no start-up transient after InitDeckAt");

    /* InitDeck at the deck positions, then held at the host positions: the start-up transient. */
    if (init_deck(&plain, g_box)) return 1;
    if (loads(plain, f_plain0, n)) return 1;
    for (int s = 1; s <= BOX_NSTEP; ++s) {
        if (step(plain, BOX_DT, q, v, a, n) || loads(plain, f, n)) return 1;
        for (int i = 12; i < n; ++i) {
            ring = fmax(ring, fabs(f[i] - f_at[i]));
            if (s == BOX_NSTEP) ring_end = fmax(ring_end, fabs(f[i] - f_at[i]));
        }
    }
    printf("InitDeck z=%.3f: line 4 fairlead load (%.4f, %.4f, %.4f) N\n", DECK_Z, f_plain0[21], f_plain0[22],
           f_plain0[23]);
    printf("InitDeck, then held at z=%.3f for 1 s: max deviation from the equilibrium %.3e N "
           "(%.3e N at 1 s)\n",
           HOST_Z, ring, ring_end);
    CHECK(ring > 1.0e4 * drift && ring > 1.0e-2, "the plain start rings");

    /* Rejections keep the previous model. */
    CableDyn_InitDeckAt(at, g_box, (int)strlen(g_box), q_host, 9, &err);
    CableDyn_GetLastError(at, msg, (int)sizeof msg);
    CHECK(err == CD_C_BAD_INPUT && strstr(msg, "coupled_positions") != NULL, "a wrong q0 length is rejected");
    CableDyn_InitDeckAt(at, g_box, (int)strlen(g_box), q_host, 13, &err);
    CHECK(err == CD_C_BAD_INPUT, "a q0 length that is not a multiple of 3 is rejected");
    CableDyn_InitDeckAt(at, g_box, (int)strlen(g_box), NULL, BOX_NHOST, &err);
    CHECK(err == CD_C_BAD_INPUT, "a NULL q0 is rejected");
    q_host[5] = NAN;
    CableDyn_InitDeckAt(at, g_box, (int)strlen(g_box), q_host, BOX_NHOST, &err);
    CHECK(err == CD_C_BAD_INPUT, "a non-finite q0 is rejected");
    CableDyn_InitDeckAt(at, g_box, (int)strlen(g_box), NULL, 0, &err);
    CHECK(err == CD_C_BAD_INPUT, "an empty q0 is rejected for a deck with Coupled points");
    CHECK(CableDyn_IsInitialized(at) && CableDyn_NCoupledDOF(at, &err) == BOX_NDOF,
          "a rejected InitDeckAt keeps the model");

    CableDyn_Close(&at, &err);
    CableDyn_Close(&full, &err);
    CableDyn_Close(&moved, &err);
    CableDyn_Close(&plain, &err);
    return 0;
}

/* A smooth host motion of the Coupled blocks: amplitude amp, sample time t. */
static void host_motion(const int *types, int nb, const double *q0, double amp, double t, double *q, double *v,
                        double *a)
{
    const double w = 2.0 * acos(-1.0);
    for (int k = 0; k < nb; ++k) {
        for (int c = 0; c < 3; ++c) {
            q[3 * k + c] = q0[3 * k + c];
            v[3 * k + c] = 0.0;
            a[3 * k + c] = 0.0;
        }
        if (types[k] != CD_C_POINT_COUPLED) continue;
        q[3 * k + 0] += amp * sin(w * t);
        v[3 * k + 0] = amp * w * cos(w * t);
        a[3 * k + 0] = -amp * w * w * sin(w * t);
        q[3 * k + 2] += 0.5 * amp * (1.0 - cos(w * t));
        v[3 * k + 2] = 0.5 * amp * w * sin(w * t);
        a[3 * k + 2] = 0.5 * amp * w * w * cos(w * t);
    }
}

/* Save, take n_trial trial steps, restore (n_restore times) and take n_steps steps; the result must
 * equal a fresh model taking the same n_steps steps, bit for bit. */
static int check_save_restore(const char *deck, const double *q_host, int n_host, double dt, double amp,
                              int n_trial, int n_steps, const char *route)
{
    void *h = NULL, *ref = NULL;
    double *q0, *q, *v, *a, *f, *f_ref, *mq, *mv, *ma, *rq, *rv, *ra;
    int err = -1, n, *ids, *types, status = 1, n_lines;
    char label[160];

    if (init_deck_at(&h, deck, q_host, n_host) || init_deck_at(&ref, deck, q_host, n_host)) return 1;
    n = CableDyn_NCoupledDOF(h, &err);
    n_lines = CableDyn_NLines(h, &err);
    q0 = calloc((size_t)n, sizeof *q0);
    q = calloc((size_t)n, sizeof *q);
    v = calloc((size_t)n, sizeof *v);
    a = calloc((size_t)n, sizeof *a);
    f = calloc((size_t)n, sizeof *f);
    f_ref = calloc((size_t)n, sizeof *f_ref);
    mq = calloc((size_t)n, sizeof *mq);
    mv = calloc((size_t)n, sizeof *mv);
    ma = calloc((size_t)n, sizeof *ma);
    rq = calloc((size_t)n, sizeof *rq);
    rv = calloc((size_t)n, sizeof *rv);
    ra = calloc((size_t)n, sizeof *ra);
    ids = calloc((size_t)(n / 3), sizeof *ids);
    types = calloc((size_t)(n / 3), sizeof *types);
    CHECK(q0 && q && v && a && f && f_ref && mq && mv && ma && rq && rv && ra && ids && types, "allocation");
    CableDyn_GetCoupledMotion(h, q0, v, a, n, &err);
    CableDyn_GetCoupledLayout(h, ids, types, n / 3, &err);
    CHECK(err == CD_C_OK, "layout");

    CableDyn_RestoreState(h, &err);
    CHECK(err == CD_C_BAD_INPUT, "restore without a save is rejected");
    CableDyn_SaveState(h, &err);
    CHECK(err == CD_C_OK, "save");
    /* Two outer correctors that take other (larger) trial motions before the final pass. */
    for (int corrector = 0; corrector < 2; ++corrector) {
        for (int s = 1; s <= n_trial; ++s) {
            host_motion(types, n / 3, q0, 2.0 * amp * (corrector + 1), s * dt, q, v, a);
            if (step(h, dt, q, v, a, n)) goto done;
        }
        CableDyn_RestoreState(h, &err);
        CHECK(err == CD_C_OK, "restore");
    }
    for (int s = 1; s <= n_steps; ++s) {
        host_motion(types, n / 3, q0, amp, s * dt, q, v, a);
        if (step(h, dt, q, v, a, n) || step(ref, dt, q, v, a, n)) goto done;
        if (loads(h, f, n) || loads(ref, f_ref, n)) goto done;
        CableDyn_GetCoupledMotion(h, mq, mv, ma, n, &err);
        CableDyn_GetCoupledMotion(ref, rq, rv, ra, n, &err);
        snprintf(label, sizeof label, "%s: restore + step %d equals the fresh model bit for bit", route, s);
        CHECK(memcmp(f, f_ref, (size_t)n * sizeof *f) == 0, label);
        CHECK(memcmp(mq, rq, (size_t)n * sizeof *mq) == 0 && memcmp(mv, rv, (size_t)n * sizeof *mv) == 0 &&
                  memcmp(ma, ra, (size_t)n * sizeof *ma) == 0,
              label);
    }
    /* Every line node agrees too. */
    for (int l = 0; l < n_lines; ++l) {
        int id = 0, nn = 0, sub = 0;
        double *xh, *xr;
        CableDyn_GetObjectInfo(h, CD_C_OBJ_LINE, l, &id, &nn, &sub, &err);
        xh = calloc((size_t)(3 * nn), sizeof *xh);
        xr = calloc((size_t)(3 * nn), sizeof *xr);
        CHECK(xh && xr, "allocation");
        CableDyn_GetLineValues(h, l, CD_C_LINE_POSITION, xh, 3 * nn, &err);
        CableDyn_GetLineValues(ref, l, CD_C_LINE_POSITION, xr, 3 * nn, &err);
        snprintf(label, sizeof label, "%s: line %d nodes after the restore equal the fresh model", route, id);
        CHECK(err == CD_C_OK && memcmp(xh, xr, (size_t)(3 * nn) * sizeof *xh) == 0, label);
        free(xh);
        free(xr);
    }
    printf("%s: %d trial steps x 2 correctors, restore, %d steps: bit-identical to the fresh model\n", route,
           n_trial, n_steps);

    /* A new initialization discards the save. */
    CableDyn_InitDeckAt(h, deck, (int)strlen(deck), q_host, n_host, &err);
    CHECK(err == CD_C_OK, "re-initialization");
    CableDyn_RestoreState(h, &err);
    CHECK(err == CD_C_BAD_INPUT, "a re-initialization discards the save");
    status = 0;

done:
    free(q0);
    free(q);
    free(v);
    free(a);
    free(f);
    free(f_ref);
    free(mq);
    free(mv);
    free(ma);
    free(rq);
    free(rv);
    free(ra);
    free(ids);
    free(types);
    CableDyn_Close(&h, &err);
    CableDyn_Close(&ref, &err);
    return status;
}

static int check_aggregate(void)
{
    void *h = NULL, *moved = NULL;
    int err = -1, n, ids[2], types[2];
    double q[MIX_NDOF], v[MIX_NDOF], a[MIX_NDOF], q_host[MIX_NDOF], f[MIX_NDOF], f_moved[MIX_NDOF];
    const char *half_moved = "capi_host_mixed_half.dat", *moved_deck = "capi_host_mixed_moved.dat";
    char msg[1025];

    if (init_deck(&h, g_mixed)) return 1;
    n = CableDyn_NCoupledDOF(h, &err);
    CHECK(err == CD_C_OK && n == MIX_NDOF, "the aggregate coupled vector holds the two fairleads only");
    CableDyn_GetCoupledLayout(h, ids, types, 2, &err);
    CHECK(err == CD_C_OK && ids[0] == 2 && ids[1] == 4, "aggregate blocks: the chain, then the cable fairlead");
    CHECK(types[0] == CD_C_POINT_COUPLED && types[1] == CD_C_POINT_COUPLED, "aggregate blocks are Coupled");
    CableDyn_GetCoupledMotion(h, q, v, a, n, &err);
    CHECK(err == CD_C_OK, "aggregate coupled motion");

    /* InitDeckAt on the aggregate: the fairleads 1 m outward and 0.5 m down. */
    memcpy(q_host, q, sizeof q);
    q_host[0] -= 1.0;
    q_host[2] -= 0.5;
    q_host[3] += 1.0;
    q_host[5] -= 0.5;
    CableDyn_InitDeckAt(h, g_mixed, (int)strlen(g_mixed), q_host, MIX_NDOF, &err);
    if (err != CD_C_OK) {
        report(h, "aggregate InitDeckAt", err);
        return 1;
    }
    CableDyn_GetCoupledMotion(h, q, v, a, n, &err);
    CHECK(err == CD_C_OK && memcmp(q, q_host, sizeof q) == 0, "aggregate fairleads start exactly at q0");
    if (loads(h, f, n)) return 1;
    printf("aggregate InitDeckAt: fairlead loads (%.1f, %.1f, %.1f) N and (%.1f, %.1f, %.1f) N\n", f[0], f[1],
           f[2], f[3], f[4], f[5]);

    /* Bit for bit the statics of the deck with both fairleads moved to q0. */
    CHECK(q_host[0] == -6.0 && q_host[2] == -5.5 && q_host[3] == 6.0 && q_host[5] == -5.5,
          "the mixed deck fairleads start at (-5, 0, -5) and (5, 0, -5)");
    if (write_moved_deck(g_mixed, half_moved, "-5.0      0.0     -5.0", "-6.0      0.0     -5.5", 1) ||
        write_moved_deck(half_moved, moved_deck, " 5.0      0.0     -5.0", " 6.0      0.0     -5.5", 1) ||
        init_deck(&moved, moved_deck) || loads(moved, f_moved, n))
        return 1;
    CHECK(memcmp(f_moved, f, sizeof f) == 0, "aggregate InitDeckAt equals InitDeck of the moved deck");
    remove(half_moved);
    remove(moved_deck);
    CableDyn_Close(&moved, &err);

    CableDyn_InitDeckAt(h, g_mixed, (int)strlen(g_mixed), q_host, 3, &err);
    CableDyn_GetLastError(h, msg, (int)sizeof msg);
    CHECK(err == CD_C_BAD_INPUT && strstr(msg, "coupled_positions") != NULL,
          "a wrong aggregate q0 length is rejected");
    CableDyn_Close(&h, &err);

    return check_save_restore(g_mixed, q_host, MIX_NDOF, MIX_DT, 0.2, 5, 10, "aggregate route");
}

int main(int argc, char **argv)
{
    double q_host[BOX_NHOST];

    if (argc != 3) {
        fputs("usage: test_c_api_host_coupling <box_points.dat> <mixed_chain_cable.dat>\n", stderr);
        return 2;
    }
    g_box = argv[1];
    g_mixed = argv[2];
    box_fairleads(HOST_Z, q_host);
    if (check_version() || check_point_layout() || check_raw_layout() || check_init_at() ||
        check_save_restore(g_box, q_host, BOX_NHOST, BOX_DT, 0.005, 7, 20, "point route") || check_aggregate()) {
        return 1;
    }
    puts("PASS: CableDyn C API host coupling");
    return 0;
}
