## Setup is similar to Vilanova et al

- Domain: `[0, 385e-6]x[0, 385e-6]x[0, 385e-6]`

- Blood visocity: `mu = 3.5e-3`

- Tissue visocity: `mu_t = 1`

- Tissue `K = 1e-15`

- `Lp = 1e-12`

- Pressue in artiry: `p1 = 50`, `p2 = 20`

- Pressue in vein: `p1 = 5`, `p2 = 10`

- Discretization of cylinder surface: `N_length = 10`, `N_angle = 10`

- Number of refinements: `5`

### Method

- To make equations well conditioned, we multiply both sides of 3-d and 1-d equation of pressure by a factor

- For 3-d, the factor is `1e+12`

- For 1-d, the factor is `1/(pi * R^3)`, where `R` is the radius of segment

We run problems with varrying radius and implicit-explicit mixed parameter `theta`.

## t_1

- Radius of artiry and vein are: `3.125e-6` and `6.25e-6`

- Mesh size of 3-d domain is `h = 385e-6 / 50`

- Aspect ratio of mesh size to minimum vessel radius is: `2.464`

- Implicit-explicit scheme parameter: `theta = 0.5`

- Result: `t_1_p.png` and `t_1_v.png`

## t_2

- Radius of artiry and vein are: `10e-6` and `20e-6`

- Rest is same as `t_1`

## t_3

- Radius of artiry and vein are: `10e-6` and `20e-6`

- Tissue visocity: `mu_t = 1.0e-3`

- Tissue conductivity: `K_t = 1.0e-18` so that ratio remains `K_t / mu_t = 1e-15`

- Rest is same as `t_1`

- Velocity calculations implemented correctly.

- As with the equations for pressure, we multiplied the both side of weak form of velocity by factor `1/Lp = 1e+12`.

- Result: `t_3_v_and_p.png`
