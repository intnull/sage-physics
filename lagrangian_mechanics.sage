def general_derivative(f, x):
    """
    Derivative of f with respect to the symbolic function x
    """
    tempX = SR.symbol()
    return f.subs_expr({x: tempX}).diff(tempX).subs_expr({tempX: x})
    
def euler_lagrange_equation(L, q):
    """
    Euler-Lagrange equation corresponding to the generalized coordinate q
    """
    return diff(general_derivative(L, diff(q, t)), t) == general_derivative(L, q)
