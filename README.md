"""

  `module Arcus`

Module Arcus provides a data type Arc, representing points on the unit circle,
or an agular value in the range `(-π,π]`.
For a point given as `(sin(α), cos(α))`, the absolutely lower value of both is stored.
Quadrant information is stored by multiplying a corresponding power of two.
For the first quadrant with `α` in  `[-π/4,π/4]`, the Float64 value of `sin(α)` is
stored unchanged.

Usage:

  `using Arcus`

  `a = Arc(radians)`

  `b = Arc(s, c)`

Trigonometric functions and arithmetic operations are available.
Multiplication and division only with `Real` values.

