import math

def dirac(x, a=0., sig=None):
    """Simple calculation of the dirac function.

    Args:
        x (float): position for function evaluation
        a (float): dirac move
        sig (float): dirac spread (default 'None' returns 0 if x != a else 1)

    Returns:
        float: value of dirac function evaluated at x
    """

    if sig == None:
        if x != a:
            result = 0
        else:
            result = 1
    else:
        result = 1 / (sig*math.pi**0.5) * math.exp(-(x-a)**2 / sig**2)

    return result

