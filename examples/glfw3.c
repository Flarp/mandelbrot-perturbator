// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <mpfr.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "perturbator.h"

static inline int max(int a, int b) {
  return a > b ? a : b;
}

static int envi(const char *name, int def) {
  const char *e = getenv(name);
  if (e) {
    return atoi(e);
  } else {
    return def;
  }
}

static double envd(const char *name, double def) {
  const char *e = getenv(name);
  if (e) {
    return atof(e);
  } else {
    return def;
  }
}

static int envr(mpfr_t out, const char *name, const char *def) {
  const char *e = getenv(name);
  if (e) {
    return mpfr_set_str(out, e,   10, MPFR_RNDN);
  } else {
    return mpfr_set_str(out, def, 10, MPFR_RNDN);
  }
}

static const char *simple_vert =
  "#version 150 core\n"
  "out vec2 texCoord;\n"
  "const vec2 t[4] = vec2[4](vec2(0.0, 0.0), vec2(1.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 1.0));\n"
  "void main() {\n"
  "  gl_Position = vec4(t[gl_VertexID] * 2.0 - 1.0, 0.0, 1.0);\n"
  "  texCoord = t[gl_VertexID] * vec2(1.0, -1.0) + vec2(0.0, 1.0);\n"
  "}\n"
  ;

static const char *simple_frag =
  "#version 150 core\n"
  "uniform sampler2D tex;\n"
  "uniform bool show_lines;\n"
  "uniform float weight;\n"
  "in vec2 texCoord;\n"
  "out vec4 colour;\n"
  "int px(vec2 tc) {\n"
  "  vec4 t = texture(tex, tc);\n"
  "  int k = 0;\n"
  "  k += (t.x / 2.0) - floor(t.x / 2.0) >= 0.5 ? 1 : 2;\n"
  "  k += t.z >= 0.0 ? 4 : 8;\n"
  "  k += 16 * int(floor(t.x));\n"
  "  return k;\n"
  "}\n"
/*
  "float sqr(float l) {\n"
  "  return l * l;\n"
  "}\n"
  "float sqr(vec2 l) {\n"
  "  return sqr(l.x + l.y);\n"
  "}\n"
*/
  "void main() {\n"
  "  float e = 1.0;\n"
  "  vec2 dx = dFdx(texCoord);\n"
  "  vec2 dy = dFdy(texCoord);\n"
  "  if (show_lines) {\n"
  "    e = px(texCoord) == px(texCoord + dx + dy) && px(texCoord + dx) == px(texCoord + dy) ? 1.0 : 0.0;\n"
  "  }\n"
  "  vec2 me = texture(tex, texCoord).xy;\n"
/*
  "  vec2 l1 = texture(tex, texCoord - dx).xy - me;\n"
  "  vec2 l2 = texture(tex, texCoord + dx).xy - me;\n"
  "  vec2 l3 = texture(tex, texCoord - dy).xy - me;\n"
  "  vec2 l4 = texture(tex, texCoord + dy).xy - me;\n"
  "  float s = weight + sqr(l1) + sqr(l2) + sqr(l3) + sqr(l4);\n"
  "  s = sqrt(weight / s);\n"
*/
  "  float s = tanh(clamp(texture(tex, texCoord).w, 0.0, 8.0));\n"
  "  colour = vec4(dot(me, me) <= 0.0 ? vec3(1.0, 0.7, 0.0) : vec3(mix(0.0, mix(0.9, 1.0, e), s)), 1.0);\n"
  "}\n"
  ;

struct state_s {
  GLFWwindow *window;

  bool should_restart;
  bool should_redraw;
  bool should_save_now;
  bool should_save_when_done;
  bool should_view_morph;

  bool show_lines;
  GLuint show_lines_u;
  double log2weight;
  GLuint weight_u;

  int width;
  int height;

  int precision;
  mpfr_t centerx;
  mpfr_t centery;
  mpfr_t radius;

};
typedef struct state_s state_t;

static void button_handler(GLFWwindow *window, int button, int action, int mods) {
  if (action == GLFW_PRESS) {
    double x = 0, y = 0;
    state_t *state = glfwGetWindowUserPointer(window);
    glfwGetCursorPos(window, &x, &y);
    mpfr_t cx, cy, r;
    mpfr_init2(cx, state->precision);
    mpfr_init2(cy, state->precision);
    mpfr_init2(r, 53);
    mpfr_set(cx, state->centerx, MPFR_RNDN);
    mpfr_set(cy, state->centery, MPFR_RNDN);
    mpfr_set(r, state->radius, MPFR_RNDN);
    double w = state->width;
    double h = state->height;
    double dx = 2 * ((x + 0.5) / w - 0.5) * (w / h);
    double dy = 2 * (0.5 - (y + 0.5) / h);
    mpfr_t ddx, ddy;
    mpfr_init2(ddx, 53);
    mpfr_init2(ddy, 53);
    mpfr_mul_d(ddx, r, dx, MPFR_RNDN);
    mpfr_mul_d(ddy, r, dy, MPFR_RNDN);
    switch (button) {
      case GLFW_MOUSE_BUTTON_LEFT:
        mpfr_mul_2si(ddx, ddx, -1, MPFR_RNDN);
        mpfr_mul_2si(ddy, ddy, -1, MPFR_RNDN);
        mpfr_add(cx, cx, ddx, MPFR_RNDN);
        mpfr_add(cy, cy, ddy, MPFR_RNDN);
        mpfr_mul_2si(r, r, -1, MPFR_RNDN);
        state->precision += 1;
        mpfr_set_prec(state->centerx, state->precision);
        mpfr_set_prec(state->centery, state->precision);
        mpfr_set(state->centerx, cx, MPFR_RNDN);
        mpfr_set(state->centery, cy, MPFR_RNDN);
        mpfr_set(state->radius, r, MPFR_RNDN);
        state->should_restart = true;
        break;
      case GLFW_MOUSE_BUTTON_RIGHT:
        mpfr_sub(cx, cx, ddx, MPFR_RNDN);
        mpfr_sub(cy, cy, ddy, MPFR_RNDN);
        mpfr_mul_2si(r, r, 1, MPFR_RNDN);
        state->precision -= 1;
        mpfr_set_prec(state->centerx, state->precision);
        mpfr_set_prec(state->centery, state->precision);
        mpfr_set(state->centerx, cx, MPFR_RNDN);
        mpfr_set(state->centery, cy, MPFR_RNDN);
        mpfr_set(state->radius, r, MPFR_RNDN);
        state->should_restart = true;
        break;
      case GLFW_MOUSE_BUTTON_MIDDLE:
        mpfr_add(cx, cx, ddx, MPFR_RNDN);
        mpfr_add(cy, cy, ddy, MPFR_RNDN);
        mpfr_set(state->centerx, cx, MPFR_RNDN);
        mpfr_set(state->centery, cy, MPFR_RNDN);
        state->should_restart = true;
        break;
    }
    mpfr_clear(cx);
    mpfr_clear(cy);
    mpfr_clear(r);
    mpfr_clear(ddx);
    mpfr_clear(ddy);
  }
(void) mods;
}

static void key_handler(GLFWwindow *window, int key, int scancode, int action, int mods) {
  state_t *state = glfwGetWindowUserPointer(window);
  if (action == GLFW_PRESS) {
    switch (key) {
      case GLFW_KEY_Q:
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;
      case GLFW_KEY_M:
        state->should_view_morph = true;
        break;
      case GLFW_KEY_L:
        state->show_lines = ! state->show_lines;
        state->should_redraw = true;
        break;
      case GLFW_KEY_S:
        if (mods & GLFW_MOD_SHIFT) {
          state->should_save_now = true;
        } else {
          state->should_save_when_done = true;
        }
        break;
      case GLFW_KEY_9:
        state->log2weight -= 0.25;
        state->should_redraw = true;
        break;
      case GLFW_KEY_0:
        state->log2weight += 0.25;
        state->should_redraw = true;
        break;
    }
  }
(void) scancode;
}

static void refresh_callback(void *user_pointer) {
  state_t *state = user_pointer;
  assert(state);
  glUniform1i(state->show_lines_u, state->show_lines);
  glUniform1f(state->weight_u, exp2f(state->log2weight));
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  glfwSwapBuffers(state->window);
}

static void handle_view_morph(struct perturbator *context, state_t *state) {
  if (state->should_view_morph) {
    state->should_view_morph = false;
    mpfr_t x, y, r;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);
    perturbator_get_primary_reference(context, x, y);
    mpfr_init2(r, 53);
    mpfr_set(r, state->radius, MPFR_RNDN);
    int e = mpfr_get_exp(r);
    mpfr_mul_2si(r, r, e / 2, MPFR_RNDN);
    int p = max(53, 53 - mpfr_get_exp(r));
    state->precision = p;
    mpfr_set_prec(state->centerx, p);
    mpfr_set_prec(state->centery, p);
    mpfr_set(state->centerx, x, MPFR_RNDN);
    mpfr_set(state->centery, y, MPFR_RNDN);
    mpfr_set(state->radius, r, MPFR_RNDN);
    state->should_restart = true;
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(r);
  }
}

extern int main(int argc, char **argv) {
  state_t state;
  memset(&state, 0, sizeof(state));

  int workers = envi("threads", 4);
  int width = envi("width", 1280);
  int height = envi("height", 720);
  int maxiters = envi("maxiters", 1 << 18);
  int chunk = envi("chunk", 1 << 8);
  double escape_radius = envd("escaperadius", 25);
  double glitch_threshold = envd("glitchthreshold", 1e-6);
  int precision = envi("precision", 53);
  int zoom = envi("zoom", 0);

  mpfr_init2(state.radius, 53);
  envr(state.radius, "radius", "2.0");

  int e = max(53, 53 - mpfr_get_exp(state.radius));
  if (e > precision) {
    fprintf(stderr, "WARNING: increasing precision to %d\n", e);
    precision = e;
  }

  mpfr_init2(state.centerx, precision);
  mpfr_init2(state.centery, precision);
  envr(state.centerx, "real", "-0.75");
  envr(state.centery, "imag", "0.0");

  glfwInit();
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
  GLFWwindow *window = glfwCreateWindow(width, height, "perturbator", 0, 0);
  glfwMakeContextCurrent(window);
  glewExperimental = GL_TRUE;
  glewInit();
  glGetError(); // discard common error from glew

  GLuint tex;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  GLuint program = glCreateProgram();
  {
    GLuint shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(shader, 1, &simple_vert, 0);
    glCompileShader(shader);
    glAttachShader(program, shader);
    glDeleteShader(shader);
  }
  {
    GLuint shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(shader, 1, &simple_frag, 0);
    glCompileShader(shader);
    glAttachShader(program, shader);
    glDeleteShader(shader);
  }
  glLinkProgram(program);
  {
    char buf[10000];
    int len = 0;
    glGetProgramInfoLog(program, 10000, &len, buf);
    if (len > 0) {
      buf[len] = 0;
      fprintf(stderr, "%s\n", buf);
    }
  }
  glUseProgram(program);
  state.show_lines_u = glGetUniformLocation(program, "show_lines");
  state.weight_u = glGetUniformLocation(program, "weight");

  uint8_t *ppm = malloc(width * height * 3);
  assert(ppm);

  state.window = window;
  state.should_restart = false;
  state.should_redraw = false;
  state.should_save_now = false;
  state.should_save_when_done = false;
  state.width = width;
  state.height = height;
  state.precision = precision;

  struct perturbator *context = perturbator_new(workers, width, height, maxiters, chunk, escape_radius, glitch_threshold);

  glfwSetWindowUserPointer(window, &state);

  if (zoom) {
    mpfr_set_d(state.radius, 256, MPFR_RNDN);
    for (int z = 0; z < zoom; ++z) {
      fprintf(stderr, "%8d FRAME\n", z);
      glfwPollEvents();
      if (glfwWindowShouldClose(window)) {
        break;
      }
      perturbator_start(context, state.centerx, state.centery, state.radius);
      perturbator_stop(context, false);
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, perturbator_get_output(context));
      refresh_callback(&state);
      glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, ppm);
      printf("P6\n%d %d\n255\n", width, height);
      for (int y = height - 1; y >= 0; --y) {
        fwrite(ppm + y * width * 3, width * 3, 1, stdout);
      }
      fflush(stdout);
      mpfr_div_2exp(state.radius, state.radius, 1, MPFR_RNDN);
    }

  } else {
    glfwSetMouseButtonCallback(window, button_handler);
    glfwSetKeyCallback(window, key_handler);

    bool first = true;
    while (! glfwWindowShouldClose(window)) {

      state.should_restart = false;
  //    fprintf(stderr, "start\n");
      if (! first) {
        perturbator_stop(context, true);
      }
      first = false;
      perturbator_start(context, state.centerx, state.centery, state.radius);

  //    fprintf(stderr, "wait_timeout\n");
      while (perturbator_active(context)) {
        struct timespec delta = { 0, 33333333 };
        nanosleep(&delta, 0);
  //      fprintf(stderr, "refresh\n");
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, perturbator_get_output(context));
        refresh_callback(&state);

        glfwPollEvents();
        if (glfwWindowShouldClose(window)) {
          break;
        }

        if (state.should_save_now) {
          state.should_save_now = false;
          glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, ppm);
          printf("P6\n%d %d\n255\n", width, height);
          fwrite(ppm, width * height * 3, 1, stdout);
          fflush(stdout);
        }

        handle_view_morph(context, &state);

        if (state.should_restart) {
          state.should_restart = false;
  //        fprintf(stderr, "start\n");
          perturbator_stop(context, true);
          perturbator_start(context, state.centerx, state.centery, state.radius);
        }
  //      fprintf(stderr, "wait_timeout\n");
      }

      if (glfwWindowShouldClose(window)) {
        break;
      }

      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, perturbator_get_output(context));
      refresh_callback(&state);

      while (! state.should_restart) {

        if (state.should_save_now || state.should_save_when_done) {
          state.should_save_now = false;
          state.should_save_when_done = false;
          glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, ppm);
          printf("P6\n%d %d\n255\n", width, height);
          fwrite(ppm, width * height * 3, 1, stdout);
          fflush(stdout);
        }

        glfwWaitEvents();
        if (glfwWindowShouldClose(window)) {
          break;
        }

  //      if (state.should_redraw) {
  //        state.should_redraw = false;
  //        fprintf(stderr, "refresh\n");
          refresh_callback(&state);
  //      }

        handle_view_morph(context, &state);

      }

    }

  //  fprintf(stderr, "delete\n");
    perturbator_stop(context, true);
  }

//  image_delete(context);
  free(ppm);
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
(void) argc;
(void) argv;
}
