load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

def bitpit_library(name, **kwargs):
    cc_library(
        name = name,
        srcs = native.glob(["src/{}/*.cpp".format(name)]),
        hdrs = native.glob(["src/{}/*.hpp".format(name)]) + ["external/LAPACKE/include/bitpit_private_lapacke.hpp"],
        includes = ["src/{}".format(name), "external/LAPACKE/include"],
        textual_hdrs = native.glob(["src/{}/*.tpp".format(name)]),
        copts = ["-Isrc/{}".format(name)],
        alwayslink = 1,
        **kwargs
    )

def bitpit_binary(name, **kwargs):
    cc_binary(
        name = name,
        **kwargs
    )
