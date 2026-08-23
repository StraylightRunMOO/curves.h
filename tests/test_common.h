#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <stdio.h>
#include <math.h>
#include <string.h>

static int g_pass = 0;
static int g_fail = 0;

#define ASSERT_TRUE(cond, msg) do { \
    if (cond) { g_pass++; } \
    else { g_fail++; fprintf(stderr, "  FAIL [%s:%d] %s\n", __func__, __LINE__, msg); } \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) \
    ASSERT_TRUE(fabsf((float)(a) - (float)(b)) < (tol), msg)

#define ASSERT_INT_EQ(a, b, msg) \
    ASSERT_TRUE((a) == (b), msg)

#define RUN_TEST(fn) do { \
    fprintf(stderr, "  %s ...\n", #fn); \
    fn(); \
} while(0)

#define TEST_SUMMARY() do { \
    fprintf(stderr, "  %d passed, %d failed\n", g_pass, g_fail); \
    return g_fail > 0 ? 1 : 0; \
} while(0)

#endif
