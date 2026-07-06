# DeconvoliSTa -- SpatialDDLS 1.0.3 <-> modern keras (tf-keras 2.16) S4 compatibility shim.
#
# SpatialDDLS's `KerasOrList` class union (used by the DeconvDLModel @model slot) only lists the
# pre-2.13 keras class name `keras.engine.sequential.Sequential`. tf-keras >= 2.13 names its model
# classes `keras.src.engine.*`, so the trained model no longer satisfies the slot and storing it fails
# with: "assignment of an object of class keras.src.engine.sequential.Sequential is not valid for
# @model ... is(value, KerasOrList) is not TRUE".
#
# No TF/Keras version reconciles this on its own: the self-contained CUDA wheels (tensorflow[and-cuda])
# only exist from TF 2.14, whose Keras already uses the keras.src.* names. So we widen the union here.
# This hook fires when SpatialDDLS attaches (library(SpatialDDLS)); it uses static class names, so no
# Python binding is required. Sourced for every `Rscript` run because it lives in Rprofile.site.
local({
  setHook(packageEvent("SpatialDDLS", "attach"), function(...) {
    cls <- c("keras.src.engine.sequential.Sequential", "keras.src.engine.functional.Functional")
    for (c in cls) if (!methods::isClass(c)) methods::setClass(c, methods::representation("VIRTUAL"))
    suppressWarnings(methods::setClassUnion("KerasOrList", c("list", cls)))
  })
})
