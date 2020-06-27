'''
  Module containing quadrature approx rules
'''
import math

import numpy as np

NUMBER_OF_NODES_GH = (1, 2, 3, 4, 5)


def create_mapping_of_nodes_and_weights(number_of_nodes):
    from scipy.special import roots_hermite
    nodes, weights = roots_hermite(number_of_nodes+1)
    return {'nodes': nodes, 'weights': weights}

def GHf(f,mu, sigma, n=5): #GHf: Gauss-Hermite quadrature for f
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
    weightsnnodes = {k: create_mapping_of_nodes_and_weights(k) 
                     for k in NUMBER_OF_NODES_GH}
    res = 1/(np.sqrt(np.pi))
    sum_ = 0
    for i in range(n+1):
        node = weightsnnodes[n]['nodes'][i]
        part_sum = mu + np.sqrt(2*(sigma**2))*node
        weight = weightsnnodes[n]['weights'][i]
        sum_ += weight*f(part_sum)
    
    sum_res = res*sum_

    return sum_res


def GLf(f,a,b,n): #GLf: Gauss-Legendre quadrature for f
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
    result = (b - a) / 2
    sum_res = 0
    weightsnnodes = {
      0 : {'weights': [2.0],'nodes':[0.0]},
      1 : {'weights': [1.0,1.0],'nodes':[-math.sqrt(1/3),math.sqrt(1/3)]},
      2 : {'weights': [5/9,8/9,5/9],'nodes':[-math.sqrt(3/5),0,math.sqrt(3/5)]},
      3 : {'weights': [0.347855, 0.652145, 0.652145, 0.347855],'nodes':[-0.861136,-0.339981,0.339981,0.861136]},
      4 : {'weights': [0.236927, 0.478629, 0.568889, 0.478629, 0.236927],'nodes':[-0.90618, -0.538469, 0, 0.538469, 0.90618]}
      #5 : {'weights': [0.171324,0.360762,0.467914,0.467914,0.360762,0.171324],'nodes':[-0.932470,-0.661209,-0.238619,0.238619,0.661209,0.932470]}
    }

    for i in range(0,len(weightsnnodes[n]['weights'])):
        part_sum = 1/2*(((b-a)*weightsnnodes[n]['nodes'][i]) + a + b)
        sum_res += (weightsnnodes[n]['weights'][i]) * (f(part_sum))

    return result*sum_res


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
    h_hat = (b - a) / n
    sum_res = 0
    for i in range(1,n-1):
        x_i = a + (i * h_hat)
        sum_res += f(x_i)
    
    x_0 = a # + (0 * h_hat) -> but this is cero
    x_n = a + (n * h_hat)
    return (h_hat/2)*(f(x_0)+f(x_n)+(2*sum_res))


# From class
def Rf(f,a,b):
    node=(a+b)/2 #middle point
    h=b-a
    return h*f(node) #polynomial of zero degree

def Rcf(f,a,b,n): #Rcf: composite rectangle method
    """
    Compute numerical approximation using rectangle or mid-point method in 
    an interval.
    Nodes are generated via formula: x_i = a+(i+1/2)h_hat for i=0,1,...,n-1 and h_hat=(b-a)/n
    Args:
        f (function): function expression of integrand
        a (float): left point of interval
        b (float): right point of interval
        n (float): number of subintervals
    Returns:
        sum_res (float): numerical approximation to integral of f in the interval a,b
    """
    h_hat=(b-a)/n
    sum_res = 0
    for i in np.arange(n):
        x = a + (i+1/2)*h_hat
        sum_res+= f(x)

    return h_hat*sum_res 


if __name__ == "__main__":
    sigma = 0.25
    mu = 0.15
    f1 = lambda t: t**2
    n = 2
    GH1_approx = GHf(f1,mu,sigma, n)
    GH1_obj = 0.15
    print(GH1_approx)