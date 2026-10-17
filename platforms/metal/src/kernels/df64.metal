/**
 * This file contains functions and classes for emulating double precision using two single precision values.  The
 * following functions are based on https://andrewthall.org/papers/df64_qf128.pdf.
 */

namespace DoubleFloat {
    inline float2 quickTwoSum(float a, float b) {
        float s = a+b;
        float e = b-(s-a);
        return float2(s, e);
    }
    inline float2 twoSum(float a, float b) {
        float s = a+b;
        float v = s-a;
        float e = (a-(s-v))+(b-v);
        return float2(s, e);
    }
    inline float2 df64add(float2 a, float2 b) {
        float2 s = twoSum(a.x, b.x);
        float2 t = twoSum(a.y, b.y);
        s.y += t.x;
        s = quickTwoSum(s.x, s.y);
        s.y += t.y;
        return quickTwoSum(s.x, s.y);
    }
    inline float2 df64diff(float2 a, float2 b) {
        float2 s = twoSum(a.x, -b.x);
        float2 t = twoSum(a.y, -b.y);
        s.y += t.x;
        s = quickTwoSum(s.x, s.y);
        s.y += t.y;
        return quickTwoSum(s.x, s.y);
    }
    inline float2 twoProd(float a, float b) {
        float x = a*b;
        float y = metal::fma(a, b, -x);
        return float2(x, y);
    }
    inline float2 df64mult(float2 a, float2 b) {
        float2 p = twoProd(a.x, b.x);
        p.y += a.x*b.y;
        p.y += a.y*b.x;
        return quickTwoSum(p.x, p.y);
    }
    inline float2 df64sqr(float2 a) {
        float2 p = twoProd(a.x, a.x);
        p.y += 2*a.x*a.y;
        return quickTwoSum(p.x, p.y);
    }
    inline float2 df64div(float2 b, float2 a) {
        float xn = 1.0f/a.x;
        float yn = b.x*xn;
        float diff = (df64diff(b, df64mult(a, float2(yn, 0)))).x;
        float2 prod = twoProd(xn, diff);
        return df64add(float2(yn, 0), prod);
    }
    inline float2 df64sqrt(float2 a) {
        float xn = metal::rsqrt(a.x);
        float yn = a.x*xn;
        float diff = (df64diff(a, float2(yn*yn, 0))).x;
        float2 prod = twoProd(xn, diff) / 2;
        return df64add(float2(yn, 0), prod);
    }
    inline bool df64eq(float2 a, float2 b) {
        return a.x == b.x && a.y == b.y;
    }
    inline bool df64neq(float2 a, float2 b) {
        return a.x != b.x || a.y != b.y;
    }
    inline bool df64lt(float2 a, float2 b) {
        return (a.x < b.x || (a.x == b.x && a.y < b.y));
    }
    inline float2 df64exp(float2 a) {
        const float thresh = 1.0e-20f*metal::exp(a.x);
        float2 s = df64add(float2(1.0f, 0.0f), a);
        float2 p = df64sqr(a);
        float m = 2.0f;
        float2 f = float2(2.0f, 0.0f);
        float2 t = p/2.0f;
        while (metal::abs(t.x) > thresh) {
            s = df64add(s, t);
            p = df64mult(p, a);
            m += 1.0f;
            f = df64mult(f, float2(m, 0.0f));
            t = df64div(p, f);
        }
        return df64add(s, t);
    }
    inline float2 df64log(float2 a) {
        float2 xi = float2(0.0f, 0.0f);
        if (!df64eq(a, float2(1.0f, 0.0f))) {
            if (a.x <= 0.0f)
                xi = float2(metal::log(a.x));
            else {
                xi.x = metal::log(a.x);
                xi = df64add(df64add(xi, df64mult(df64exp(-xi), a)), float2(-1.0f, 0.0f));
            }
        }
        return xi;
    }
    inline float2 df64sin(float2 a) {
        const float thresh = 1.0e-20f*metal::abs(a.x);
        if (a.x == 0.0f)
            return float2(0.0f, 0.0f);
        float2 x = -df64sqr(a);
        float2 s = a;
        float2 p = a;
        float m = 1.0f;
        float2 f = float2(1.0f, 0.0f);
        while (true) {
            p = df64mult(p, x);
            m += 2.0f;
            f = df64mult(f, float2(m*(m-1), 0.0f));
            float2 t = df64div(p, f);
            s = df64add(s, t);
            if (metal::abs(t.x) < thresh)
                break;
        }
        return s;
    }
    inline float2 df64cos(float2 a) {
        const float thresh = 1.0e-20f*metal::abs(a.x);
        const float2 one = float2(1.0f, 0.0f);
        if (a.x == 0.0f)
            return one;
        float2 x = -df64sqr(a);
        float2 c = one;
        float2 p = one;
        float m = 0.0f;
        float2 f = one;
        while (true) {
            p = df64mult(p, x);
            m += 2.0f;
            f = df64mult(f, float2(m*(m-1), 0.0f));
            float2 t = df64div(p, f);
            c = df64add(c, t);
            if (metal::abs(t.x) < thresh)
                break;
        }
        return c;
    }
}

/**
 * This class represents a double-float value.  It attempts to act as a drop in replacement for double.  Note, however,
 * that while its precision is similar to a double, its range is only that of a float.
 */
class df64 {
public:
    float2 v;
    df64() = default;
    df64() threadgroup = default;
    df64(float f) {
        v = float2(f, 0.0f);
    }
    df64(int i) {
        v = float2(i, 0.0f);
    }
    df64(long i) {
        v = float2(i, 0.0f);
    }
    df64(float2 f) {
        v = f;
    }
    df64 operator-() const {
        return df64(float2(-v.x, -v.y));
    }
#define DEFINE_OPERATORS(AS) \
    template <typename T, typename = metal::enable_if_t<metal::is_same<T, float>::value>> \
    operator T() const AS { \
        return v.x; \
    } \
    AS df64& operator+=(const df64 rhs) AS { \
        v = DoubleFloat::df64add(v, rhs.v); \
        return *this; \
    } \
    AS df64& operator-=(const df64 rhs) AS { \
        v = DoubleFloat::df64diff(v, rhs.v); \
        return *this; \
    } \
    AS df64& operator*=(const df64 rhs) AS { \
        v = DoubleFloat::df64mult(v, rhs.v); \
        return *this; \
    } \
    AS df64& operator/=(const df64 rhs) AS { \
        v = DoubleFloat::df64div(v, rhs.v); \
        return *this; \
    }
    DEFINE_OPERATORS(thread)
    DEFINE_OPERATORS(threadgroup)
    DEFINE_OPERATORS(device)
#undef DEFINE_OPERATORS
};

/**
 * Define operators that act on two df64 values.
 */
df64 operator+(df64 a, df64 b) {
    return df64(DoubleFloat::df64add(a.v, b.v));
}
df64 operator-(df64 a, df64 b) {
    return df64(DoubleFloat::df64diff(a.v, b.v));
}
df64 operator*(df64 a, df64 b) {
    return df64(DoubleFloat::df64mult(a.v, b.v));
}
df64 operator/(df64 a, df64 b) {
    return df64(DoubleFloat::df64div(a.v, b.v));
}
bool operator==(df64 a, df64 b) {
    return DoubleFloat::df64eq(a.v, b.v);
}
bool operator!=(df64 a, df64 b) {
    return DoubleFloat::df64neq(a.v, b.v);
}
bool operator<(df64 a, df64 b) {
    return DoubleFloat::df64lt(a.v, b.v);
}
bool operator>(df64 a, df64 b) {
    return b < a;
}
bool operator<=(df64 a, df64 b) {
    return !(a > b);
}
bool operator>=(df64 a, df64 b) {
    return !(a < b);
}

/**
 * Define comparison operators that compare a df64 to another value of a different type.
 */
#define DEFINE_MIXED_COMPARISONS(OP) \
template <typename T, typename = metal::enable_if_t<metal::is_arithmetic<T>::value>> \
inline bool operator OP(df64 a, T b) { \
    return a OP df64(b); \
} \
template <typename T, typename = metal::enable_if_t<metal::is_arithmetic<T>::value>> \
inline bool operator OP(T a, df64 b) { \
    return df64(a) OP b; \
}
DEFINE_MIXED_COMPARISONS(==)
DEFINE_MIXED_COMPARISONS(!=)
DEFINE_MIXED_COMPARISONS(<)
DEFINE_MIXED_COMPARISONS(>)
DEFINE_MIXED_COMPARISONS(<=)
DEFINE_MIXED_COMPARISONS(>=)
#undef DEFINE_MIXED_COMPARISONS

/**
 * This is a vector type containing two df64 values.  It attempts to act as a replacement for double2.
 */
struct df642 {
    df64 x, y;
    df642(df64 a, df64 b) : x(a), y(b) {
    }
    explicit df642(df64 a) : x(a), y(a) {
    }
#define DEFINE_OPERATORS2(AS) \
    template <class T> \
    AS df642& operator+=(T rhs) AS { \
        df64 d(rhs); \
        x += d; \
        y += d; \
        return *this; \
    } \
    template <> \
    AS df642& operator+=(df642 rhs) AS { \
        x += rhs.x; \
        y += rhs.y; \
        return *this; \
    } \
    template <class T> \
    AS df642& operator-=(T rhs) AS { \
        df64 d(rhs); \
        x -= d; \
        y -= d; \
        return *this; \
    } \
    template <> \
    AS df642& operator-=(df642 rhs) AS { \
        x -= rhs.x; \
        y -= rhs.y; \
        return *this; \
    } \
    template <class T> \
    AS df642& operator*=(T rhs) AS { \
        df64 d(rhs); \
        x *= d; \
        y *= d; \
        return *this; \
    } \
    template <> \
    AS df642& operator*=(df642 rhs) AS { \
        x *= rhs.x; \
        y *= rhs.y; \
        return *this; \
    } \
    template <class T> \
    AS df642& operator/=(T rhs) AS { \
        df64 d(rhs); \
        x /= d; \
        y /= d; \
        return *this; \
    } \
    template <> \
    AS df642& operator/=(df642 rhs) AS { \
        x /= rhs.x; \
        y /= rhs.y; \
        return *this; \
    }
    DEFINE_OPERATORS2(thread)
    DEFINE_OPERATORS2(threadgroup)
    DEFINE_OPERATORS2(device)
#undef DEFINE_OPERATORS2
};

/**
 * This is a vector type containing three df64 values.  It attempts to act as a replacement for double3.
 */
struct df643 {
    df64 x, y, z;
    df643(df64 a, df64 b, df64 c) : x(a), y(b), z(c) {
    }
    explicit df643(df64 a) : x(a), y(a), z(a) {
    }
#define DEFINE_OPERATORS3(AS) \
    template <class T> \
    AS df643& operator+=(T rhs) AS { \
        df64 d(rhs); \
        x += d; \
        y += d; \
        z += d; \
        return *this; \
    } \
    template <> \
    AS df643& operator+=(df643 rhs) AS { \
        x += rhs.x; \
        y += rhs.y; \
        z += rhs.z; \
        return *this; \
    } \
    template <class T> \
    AS df643& operator-=(T rhs) AS { \
        df64 d(rhs); \
        x -= d; \
        y -= d; \
        z -= d; \
        return *this; \
    } \
    template <> \
    AS df643& operator-=(df643 rhs) AS { \
        x -= rhs.x; \
        y -= rhs.y; \
        z -= rhs.z; \
        return *this; \
    } \
    template <class T> \
    AS df643& operator*=(T rhs) AS { \
        df64 d(rhs); \
        x *= d; \
        y *= d; \
        z *= d; \
        return *this; \
    } \
    template <> \
    AS df643& operator*=(df643 rhs) AS { \
        x *= rhs.x; \
        y *= rhs.y; \
        z *= rhs.z; \
        return *this; \
    } \
    template <class T> \
    AS df643& operator/=(T rhs) AS { \
        df64 d(rhs); \
        x /= d; \
        y /= d; \
        z /= d; \
        return *this; \
    } \
    template <> \
    AS df643& operator/=(df643 rhs) AS { \
        x /= rhs.x; \
        y /= rhs.y; \
        z /= rhs.z; \
        return *this; \
    }
    DEFINE_OPERATORS3(thread)
    DEFINE_OPERATORS3(threadgroup)
    DEFINE_OPERATORS3(device)
#undef DEFINE_OPERATORS3
};

/**
 * This is a vector type containing four df64 values.  It attempts to act as a replacement for double4.
 */
struct df644 {
    df64 x, y, z, w;
    df644(df64 a, df64 b, df64 c, df64 d) : x(a), y(b), z(c), w(d) {
    }
    explicit df644(df64 a) : x(a), y(a), z(a), w(a) {
    }
#define DEFINE_OPERATORS4(AS) \
    template <class T> \
    AS df644& operator+=(T rhs) AS { \
        df64 d(rhs); \
        x += d; \
        y += d; \
        z += d; \
        w += d; \
        return *this; \
    } \
    template <> \
    AS df644& operator+=(df644 rhs) AS { \
        x += rhs.x; \
        y += rhs.y; \
        z += rhs.z; \
        w += rhs.w; \
        return *this; \
    } \
    template <class T> \
    AS df644& operator-=(T rhs) AS { \
        df64 d(rhs); \
        x -= d; \
        y -= d; \
        z -= d; \
        w -= d; \
        return *this; \
    } \
    template <> \
    AS df644& operator-=(df644 rhs) AS { \
        x -= rhs.x; \
        y -= rhs.y; \
        z -= rhs.z; \
        w -= rhs.w; \
        return *this; \
    } \
    template <class T> \
    AS df644& operator*=(T rhs) AS { \
        df64 d(rhs); \
        x *= d; \
        y *= d; \
        z *= d; \
        w *= d; \
        return *this; \
    } \
    template <> \
    AS df644& operator*=(df644 rhs) AS { \
        x *= rhs.x; \
        y *= rhs.y; \
        z *= rhs.z; \
        w *= rhs.w; \
        return *this; \
    } \
    template <class T> \
    AS df644& operator/=(T rhs) AS { \
        df64 d(rhs); \
        x /= d; \
        y /= d; \
        z /= d; \
        w /= d; \
        return *this; \
    } \
    template <> \
    AS df644& operator/=(df644 rhs) AS { \
        x /= rhs.x; \
        y /= rhs.y; \
        z /= rhs.z; \
        w /= rhs.w; \
        return *this; \
    }
    DEFINE_OPERATORS4(thread)
    DEFINE_OPERATORS4(threadgroup)
    DEFINE_OPERATORS4(device)
#undef DEFINE_OPERATORS4
};

/**
 * Define operators that act on two vectors of df64 values.
 */
inline df642 operator+(df642 a, df642 b) {
    return df642(a.x+b.x, a.y+b.y);
}
inline df642 operator-(df642 a, df642 b) {
    return df642(a.x-b.x, a.y-b.y);
}
inline df642 operator*(df642 a, df642 b) {
    return df642(a.x*b.x, a.y*b.y);
}
inline df642 operator/(df642 a, df642 b) {
    return df642(a.x/b.x, a.y/b.y);
}
inline df643 operator+(df643 a, df643 b) {
    return df643(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline df643 operator-(df643 a, df643 b) {
    return df643(a.x-b.x, a.y-b.y, a.z-b.z);
}
inline df643 operator*(df643 a, df643 b) {
    return df643(a.x*b.x, a.y*b.y, a.z*b.z);
}
inline df643 operator/(df643 a, df643 b) {
    return df643(a.x/b.x, a.y/b.y, a.z/b.z);
}
inline df644 operator+(df644 a, df644 b) {
    return df644(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}
inline df644 operator-(df644 a, df644 b) {
    return df644(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}
inline df644 operator*(df644 a, df644 b) {
    return df644(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
}
inline df644 operator/(df644 a, df644 b) {
    return df644(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
}

/**
 * Define operators that act on two values, one being a df64 scalar or vector and the other being a primitive type.
 */
#define DEFINE_BINARY_OPERATORS(T, U) \
inline T operator+(T a, U b) { \
    return a+T(b); \
} \
inline T operator+(U a, T b) { \
    return T(a)+b; \
} \
inline T operator-(T a, U b) { \
    return a-T(b); \
} \
inline T operator-(U a, T b) { \
    return T(a)-b; \
} \
inline T operator*(T a, U b) { \
    return a*T(b); \
}\
inline T operator*(U a, T b) { \
    return T(a)*b; \
} \
inline T operator/(T a, U b) { \
    return a/T(b); \
} \
inline T operator/(U a, T b) { \
    return T(a)/b; \
} \

DEFINE_BINARY_OPERATORS(df64, float)
DEFINE_BINARY_OPERATORS(df64, int)
DEFINE_BINARY_OPERATORS(df64, long)
DEFINE_BINARY_OPERATORS(df642, float)
DEFINE_BINARY_OPERATORS(df642, int)
DEFINE_BINARY_OPERATORS(df642, long)
DEFINE_BINARY_OPERATORS(df643, float)
DEFINE_BINARY_OPERATORS(df643, int)
DEFINE_BINARY_OPERATORS(df643, long)
DEFINE_BINARY_OPERATORS(df644, float)
DEFINE_BINARY_OPERATORS(df644, int)
DEFINE_BINARY_OPERATORS(df644, long)
#undef DEFINE_BINARY_OPERATORS

/**
 * Math functions that act on df64 values.
 */
inline df64 sqrt(df64 a) {
    return df64(DoubleFloat::df64sqrt(a.v));
}
inline df64 rsqrt(df64 a) {
    return df64(1.0f)/df64(DoubleFloat::df64sqrt(a.v));
}
inline df64 log(df64 a) {
    return df64(DoubleFloat::df64log(a.v));
}
inline df64 exp(df64 a) {
    return df64(DoubleFloat::df64exp(a.v));
}
inline df64 sin(df64 a) {
    return df64(DoubleFloat::df64sin(a.v));
}
inline df64 cos(df64 a) {
    return df64(DoubleFloat::df64cos(a.v));
}
inline df64 tan(df64 a) {
    float2 s = DoubleFloat::df64sin(a.v);
    float2 c = DoubleFloat::df64cos(a.v);
    return df64(s)/df64(c);
}
inline df64 fabs(df64 v) {
    return (v > df64(0.0f) ? v : -v);
}
inline df643 cross(df643 a, df643 b) {
    return df643(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}
inline df644 cross(df644 a, df644 b) {
    return df644(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x, df64(0.0f));
}
inline df643 trimTo3(df644 v) {
    return df643(v.x, v.y, v.z);
}

/**
 * If we're using the fast:: versions of math functions, make sure the df64 versions will still be found.
 */
namespace metal {
    namespace fast {
        inline df64 sqrt(df64 a) {
            return ::sqrt(a);
        }
        inline df64 rsqrt(df64 a) {
            return ::rsqrt(a);
        }
        inline df64 divide(float a, df64 b) {
            return a/b;
        }
        inline df64 exp(df64 a) {
            return ::exp(a);
        }
        inline df64 log(df64 a) {
            return ::log(a);
        }
    }
}