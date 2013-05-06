glmer2stan
==========

Define Stan models using glmer-style (lme4) formulas.

Status
==========

Current state: Production usable. All models pass tests, DIC computations validated. Basic prediction function added, but still needs work. WAIC is in and working for single-formula models, but still needs some more thorough numerical testing, to find corner cases.

Horizon: Code will eventually need to be refactored, to reduce redundancy and clean up some of the compilation algorithm. Speed of Stan code could be improved during that refactoring.

Examples
==========

[Add a few examples]
