# Note: this module inserts the variable t in the global namespace!
# Conventions:
# t = time
# n = degrees of freedom
# s = number of particles
# q = [q1, ..., qn] generalized coordinates
# p = [p1, ..., pn] canonical momentum conjugate to q
# r = [r1, ..., rs] position vectors
# v = [v1, ..., vs] velocity vectors
# m = [m1, ..., ms] masses

from sage.all import SR, var, function, diff, solve
t = var('t')

def dot(f):
    r"""
    The derivative of `f` with respect to time.
    """
    return diff(f, t)

def formal_derivative(f, x):
    r"""
    The formal derivative of `f` with respect to the symbolic function `x`.
    """
    tempX = SR.symbol()
    return f.subs_expr({x: tempX}).diff(tempX).subs_expr({tempX: x})
    
def dynamical_var(s):
    r"""
    Create a formal symbolic function of `t` with the name s.
    """
    G = globals()
    if ',' in s:
        L = [i.strip() for i in s.split(',')]
        for var in L:
            G[var] = function(var, t)
            var = G[var]
        return tuple(L)
    elif ' ' in s:
        L = [i.strip() for i in s.split(' ')]
        for var in L:
            G[var] = function(var, t)
            var = G[var]
        return tuple(L)
    else:
        G[s] = function(s, t)
        return G[s] 

def kinetic_energy(v, m):
    r"""
    The kinetic energy.

    EXAMPLES:
        sage: m = var('m')
        sage: q = dynamical_var('q')
        sage: kinetic_energy(dot(q), m)
        1/2*m*D[0](q)(t)^2
    """
    try:
        v[0][0] # test if v is a vector of a vector
        sum = 0
        for s in range(len(v)):
            sum += m[s]/2 * (v[s] * v[s])
        return sum
    except:
        return m/2 * (v * v)

def euler_lagrange_equation(L, q):
    r"""
    The Euler-Lagrange equation corresponding to the generalized coordinate `q`.
    """
    try:
        n = len(q)
        result = []
        for i in range(n):
            result.append(diff(formal_derivative(L, dot(q[i])), t) == formal_derivative(L, q[i]))
        return result
    except TypeError:
        return diff(formal_derivative(L, dot(q)), t) == formal_derivative(L, q)

def poisson_bracket(f, g, q, p):
    r"""
    The poisson bracket of `f` and `g`.
    """
    try:
        n = len(q)
        sum = 0
        for i in range(n):
            sum += formal_derivative(f, q[i]) * formal_derivative(g, p[i]) - \
                formal_derivative(f, p[i]) * formal_derivative(g, q[i])
        return sum
    except TypeError:
        return formal_derivative(f, q) * formal_derivative(g, p) - \
            formal_derivative(f, p) * formal_derivative(g, q)

def hamilton_equations(H, q, p):
    r"""
    The Hamilton equations.
    """
    try:
        n = len(q)
        result = []
        for i in range(n):
            result.append(dot(q[i]) == formal_derivative(H, p[i]))
        for i in range(n):
            result.append(dot(p[i]) == - formal_derivative(H, q[i]))
        return result
    except TypeError:
        return [dot(q) == formal_derivative(H, p), dot(p) == - formal_derivative(H, q)]

def legendre_transformation(L, qdot, p):
    r"""
    The Legendre transformation of the Lagrangian `L` with respect to `qdot`.
    """
    # Warning: n=1 only so far
    eqn = p == formal_derivative(L, qdot)
    new_qdot = solve(eqn, qdot)[0].rhs()
    new_L = L.substitute({qdot: new_qdot})
    H = new_qdot * p - new_L
    return H
