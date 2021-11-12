import os
import dolfin as df

N = 10
mesh = df.RectangleMesh(df.Point(0, 0), df.Point(1, 1), N, N)

PEl = df.FiniteElement('P', mesh.ufl_cell(), 2)
REl = df.FiniteElement('R', mesh.ufl_cell(), 0)
MixedEl = df.MixedElement([PEl, REl])
V = df.FunctionSpace(mesh, MixedEl)

volume = 1

p, llambda = df.TrialFunctions(V)
q, mu = df.TestFunctions(V)

p_sol_expr = df.Expression('x[0]*x[0] + 2*x[0] - 1', degree=2)
g = df.Expression(('2*x[0]+ 2', 0), degree=1)

n = df.FacetNormal(mesh)

a = 0
a += df.inner(df.grad(p), df.grad(q)) * df.dx + llambda * q * df.dx + p * q * df.dx
a += 1./volume * llambda * mu * df.dx - 1./volume * p * mu * df.dx
b = df.Expression('-2 + 1./3. + x[0]*x[0] + 2*x[0] - 1', degree=2) * q * df.dx + df.Constant(0) * mu * df.dx
b += + df.dot(g, n) * q * df.ds

u = df.Function(V, name='u')

df.solve(a == b, u, [])

output_folder = "tmp"
file = df.File(os.path.join(output_folder, "averages.pvd"), "compressed")
file << u.sub(0)
file = df.File(os.path.join(output_folder, "averages_sol.pvd"), "compressed")
p_sol_interp = df.interpolate(p_sol_expr, V.sub(0).collapse())
p_sol_interp.rename('u_sol', '')
file << p_sol_interp

print('error = {}'.format(df.errornorm(p_sol_expr, u.split(True)[0])))
print('multiplier average = {}'.format(u.split(True)[1].vector()[:]))
print('real average = {}'.format(df.assemble(u.split()[0] * df.dx) / volume))
