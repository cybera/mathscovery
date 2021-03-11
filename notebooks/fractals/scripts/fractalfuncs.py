
import numpy as np

def NewtonRaphson(f, fprime, x0, max_iter = 100, prec = 1e-10, verbose=False):
    
    x = [x0]
    x_val = x0
        
    for i in range(max_iter):
        
        val_i = f(x[i])
        df = fprime(x[i])
        epsilon = val_i/df
        x_val = x_val - epsilon
        x.append(x_val)
        
        if abs(epsilon) < prec:
            print(abs(df))
            return x_val, x
    if verbose:
        print(f"This calculation did not converge after {max_iter} iterations")
    
    return x_val, x


def NewtonRaphsonFact(f=None, fprime=None, x0=None, max_iter = 100, prec = 1e-10, verbose=False, **kwargs):
    
    x = [x0]
    x_val = x0

    for i in range(max_iter):
        
        val_i = f(x[i])
        df = fprime(x[i])
        epsilon = val_i/df
        x_val = x_val - epsilon
        x.append(x_val)
        
        if abs(epsilon) < prec:
            return i
    if verbose:
        print(f"This calculation did not converge after {max_iter} iterations")
    
    return i

def mandelbrot(c, max_iter = 80, **kwargs):
    z = 0
    n = 0
    while abs(z) <= 2 and n < max_iter:
        z = z**2 + c
        n += 1
    return n

def CreateImageMap(function, function_args, bounds, max_iter = 80, width = 200, height = 200):
    
    re_start, re_end, im_start, im_end = bounds
    if width > 1000:
        print(f'width of {width} is too large. Your computer only has so many pixels.')
        print("try zooming in with a smaller boundary to observe more detail")
        return
        
    if height > 1000:
        print(f'height of {height} is too large. Your computer only has so many pixels.')
        print("try zooming in with a smaller boundary to observe more detail")
        return
    
    X = np.zeros([width, height])
    for x in range(0, width):
        for y in range(0,height):
            real = re_start + (x/width) * ( re_end - re_start)
            imaginary = im_start + (y/width) * ( im_end - im_start)
            guess = complex(real, imaginary)
            try:
                function_args['mult']
            except KeyError:
                function_args['mult'] = 1.1
            
            function_args['x0'] = guess
            function_args['x1'] = function_args['mult'] * guess
            function_args['c'] = guess
            try:
                m = function(**function_args)
            except Exception as e:
                print(e)
                m = max_iter
            color = 255 - int(m * 255 / max_iter)
            
            X[x,y] = color
            
    return X

def secantfact(function, x0, x1, max_iter = 200, prec = 1e-5, verbose = False, **kwargs):
    
    f1 = function(x0)
    f = function(x1)
    
    if abs(f1) < abs(f):
        rts = x0
        x1 = x1
        f1, f = f, f1
    else:
        rts = x1
    for i in range(max_iter):
        dx = (x0 - rts) * f / (f - f1)
        x0 = rts
        f1 = f
        rts += dx
        f = function(rts)
        
        if abs(dx) < prec or f == 0:
            return i
   
    if verbose:
         print(f"This calculation did not converge after {max_iter} iterations")
    return i


def schroderfact(derivative, function, secondder, x0, prec = 1e-5, max_iter = 50, **kwargs):
    check = 1
    xim1 = x0
    n = 0
    while prec < check:
        num = (function(xim1) * derivative(xim1))
        dem = derivative(xim1)**2 - function(xim1) * secondder(xim1)
        xi = xim1 - num / dem
        n += 1
        check = abs(xi - xim1)
        xim1 = xi
        # print(check)
        if n == max_iter:
            break
    return  n 

def halleyfact(derivative, function,seconder, x0, prec = 1e-5, max_iter = 100, **kwargs):
    check = 1
    xim1 = x0
    n = 0
    for i in range(max_iter):
        
        a = function(xim1) * seconder(xim1)
        b = 2 * derivative(xim1) ** 2
        c = 1 - a/b

        d = derivative(xim1) * c
        xi = xim1 - function(xim1) / d
        n += 1
        check = abs(xi - xim1)
        xim1 = xi
        eps = abs(function(xim1) / d)
        
        if eps < prec:
            break

    
    return  n 

def nderiv(func, x, eps=np.sqrt(np.finfo(float).eps)):
    ''' Takes a vector of each component of a multivariate 
    function and returns its Jacobian matrix as computed by 
    finite differences'''
    
    N = len(x)

    J = [[None for i in range(N)] for j in range(N)]
    #xh is x + h in the derivative formula 
    xh = x
    # very bad quick fix, will not work as a gradient 
    if isinstance(x[0], complex):
        eps = complex(eps, eps)
    for i in range(N):
        temp = xh[i]
        h = eps * temp
        if h == 0: h = eps
        xh[i] = temp + h # scootch that point over 
        # evalueat f(x+h)
        f = func(xh)
        
        xh[i] = temp
        fvec = func(x)
        # forward difference formula
        for j in range(N):
            J[j][i] = (f[j] - fvec[j])/h

    if len(x) == 1:
        return J[0][0]
            
    return J 
