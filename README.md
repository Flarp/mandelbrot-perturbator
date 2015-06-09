perturbator
===========

efficient deep zooming for Mandelbrot sets
Copyright (C) 2015 Claude Heiland-Allen
License GPL3+ http://www.gnu.org/licenses/gpl.html


quickstart
----------

    make -C src
    make -C examples
    ./examples/perturbator-glfw3


rendering methods
-----------------

different rendering methods used for different depths:

    kind    number maxdepth   
    plain   f32    1e-7
    plain   f64    1e-15
    perturb f32    1e-38
    perturb f64    1e-308
    perturb f80    1e-9864
    perturb f??i16 1e-9864
    perturb f??i32 1e-646456993
    perturb f??i64 1e-2776511644261678080

implemented so far:

    perturb f64 "z^2 + c" with order 24 series approximation


future api
----------

    struct buffer
    buffer      *create(int width, int height)
    void         destroy(buffer *)
    const float *get_data(const buffer *)
    int          get_width(const buffer *)
    int          get_height(const buffer *)

    struct view
    view        *create()
    void         destroy(view *)
    void         set_x(view *, const mpfr_t)
    void         set_y(view *, const mpfr_t)
    void         set_r(view *, const mpfr_t)
    void         get_x(const view *, mpfr_t)
    void         get_y(const view *, mpfr_t)
    void         get_r(const view *, mpfr_t)

    struct options
    options     *create()
    void         destroy(options *)
    void         set_threads(options *, int)
    void         set_maxiters(options *, int)
    void         set_maxperiod(options *, int)
    void         set_escape_radius(options *, double)
    void         set_glitch_radius(options *, double)

    struct context
    context     *create(enum formula, int order)
    void         destroy(context *)
    void         start(context *, buffer *, const view *, const options *)
    bool         done(context *)
    void         wait(context *)
    void         stop(context *)