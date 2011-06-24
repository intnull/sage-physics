# analytical_mechanics.sage
# Conventions:
# t = time
# n = degrees of freedom
# s = number of particles
# q = [q1, ..., qn] generalized coordinates
# p = [p1, ..., pn] canonical momentum conjugate to q
# r = [r1, ..., rs] position vectors
# v = [v1, ..., vs] velocity vectors
# m = [m1, ..., ms] masses

from sage.all import SR, var, function, diff
var('t, m')

def dot(f):
    r"""
    The derivative of f with respect to time.
    """
    return diff(f, t)

def formal_derivative(f, x):
    r"""
    The formal derivative of f with respect to the symbolic function x.
    """
    tempX = SR.symbol()
    return f.subs_expr({x: tempX}).diff(tempX).subs_expr({tempX: x})
    
def dynamical_var(s):
    r"""
    Create a formal symbolic function of t with the name s.
    """
    G = globals()
    if ',' in s:
        L = [i.strip() for i in s.split(',')]
        for i in range(len(L)):
            G[L[i]] = function(L[i], t)
            L[i] = G[L[i]]
        return tuple(L)
    elif ' ' in s:
        L = [i.strip() for i in s.split(' ')]
        for i in range(len(L)):
            G[L[i]] = function(L[i], t)
            L[i] = G[L[i]]
        return tuple(L)
    else:
        G[s] = function(s, t)
        return G[s] 

def kinetic_energy(v, m=m):
    r"""
    The kinetic energy.
    """
    try:
        k = len(v)
        sum = 0
        for i in range(k):
            sum += m[i]/2 * v[i]^2
        return sum
    except TypeError:
        return m/2 * v^2

dynamical_var('q, p')
def euler_lagrange_equation(L, q=q):
    r"""
    The Euler-Lagrange equation corresponding to the generalized coordinate q.
    """
    try:
        n = len(q)
        result = []
        for i in range(n):
            result.append(diff(formal_derivative(L, diff(q[i], t)), t) == formal_derivative(L, q[i]))
        return result
    except TypeError:
        return diff(formal_derivative(L, diff(q, t)), t) == formal_derivative(L, q)

def poisson_bracket(f, g, q=q, p=p):
    r"""
    The poisson bracket of f and g.
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

def hamilton_equations(H, q=q, p=p):
    r"""
    The Hamilton equations.
    """
    try:
        n = len(q)
        result = []
        for i in range(n):
            result.append(diff(q[i], t) == formal_derivative(H, p[i]))
        for i in range(n):
            result.append(diff(p[i], t) == - formal_derivative(H, q[i]))
        return result
    except TypeError:
        return [diff(q, t) == formal_derivative(H, p), diff(p, t) == - formal_derivative(H, q)]
