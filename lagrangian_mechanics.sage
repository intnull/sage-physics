def general_derivative(f, x):
    r"""
    The derivative of f with respect to the symbolic function x.
    """
    tempX = SR.symbol()
    return f.subs_expr({x: tempX}).diff(tempX).subs_expr({tempX: x})
    
def euler_lagrange_equation(L, q):
    r"""
    The Euler-Lagrange equation corresponding to the generalized coordinate q.
    """
    return diff(general_derivative(L, diff(q, t)), t) == general_derivative(L, q)

def dynamical_var(s):
    r"""
    Create a symbolic function of t with the name s.
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
