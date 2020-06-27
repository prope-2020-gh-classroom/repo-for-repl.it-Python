import math

import numpy as np #thir party modules

#import mymodule

def Tcf(f,a,b,n): #Tcf: composite trapezoidal method for f
    """
    Compute numerical approximation using trapezoidal method in 
    an interval.
    Nodes are generated via formula: x_i = a+ih_hat for i=0,1,...,n and h_hat=(b-a)/n
    Args:
        f (function): function expression of integrand
        a (float): left point of interval
        b (float): right point of interval
        n (float): number of subintervals
    Returns:
        sum_res (float): numerical approximation to integral of f in the interval a,b
    """
    h_hat=(b-a)/n
    acc = 0
    sum_res = 0
    for i in np.arange(n):
        x = a + i*h_hat
        acc+= f(x)
    sum_res = (0.5*h_hat)*(f(a)+f(x)+2*acc)
    return sum_res


def GLf(f,a,b,n): 
    """
    Compute numerical approximation using quadrature Gauss-Legendre.
    Weights and nodes are obtained with table for n=0,1,2,3,4
    Args:
        f (function): function expression of integrand
        a (float): left point of interval
        b (float): right point of interval
        n (float): number of subintervals
    Returns:
        sum_res (float): numerical approximation to integral of f in the interval a,b
    """
    sum_res = 0
    if n == 0:
        nodes = np.array([0])
        weights = np.array([2])
    if n == 1:
        nodes = np.array([-math.sqrt(1/3),math.sqrt(1/3)])
        weights = np.array([1,1])
    if n == 2:
        nodes = np.array([-math.sqrt(3/5),0,math.sqrt(3/5)])
        weights = np.array([5/9,8/9,5/9])
    if n == 3:
        nodes = np.array([-0.861136,-0.339981,0.339981,0.861136])
        weights = np.array([0.347855,0.652145,0.652145,0.347855])
    if n == 4:
        nodes = np.array([-0.90618,-0.538469,0,0.538469,0.90618])
        weights = np.array([0.236927,0.478629,0.568889,0.478629,0.236927])
    
    constant_change_of_variable = (b-a)/2
    for weight,node in zip(weights,nodes):
      node_change_of_variable = 0.5*((b-a)*node+a+b)
      sum_res += weight*f(node_change_of_variable)
    sum_res = constant_change_of_variable*sum_res
    return sum_res


def GHf(f,mu, sigma): #GHf: Gauss-Hermite quadrature for f
    """
    Compute numerical approximation using quadrature Gauss-Hermite.
    Weights and nodes are obtained with table in Kiusalaas for n=6
    Args:
        f (function): function expression of integrand
        mu (float): mean
        sigma (float): standard deviation
    Returns:
        sum_res (float): numerical approximation to integral of f in the interval a,b
    """
    nodes = [-2.350605,-1.335849,-0.436077,0.436077,1.335849,2.350605]
    weights = [0.453001e-2,0.157067,0.724629,0.724629,0.157067,0.453001e-2]
    sum_res = 0
    n = 6
    constant_change_of_variable = 1/math.sqrt(np.pi)
    for i in np.arange(n):
      node_i = nodes[i]
      weight_i = weights[i]
      node_change_of_variable = math.sqrt(2)*sigma* node_i + mu
      sum_res +=  weight_i * f(node_change_of_variable)
    sum_res = constant_change_of_variable*sum_res
    return sum_res  

if __name__ == "__main__":
    sigma = 0.25
    mu = 0.15
    f = lambda t: t
    GH_approx = GHf(f, mu, sigma)
    #print(GH_approx)

