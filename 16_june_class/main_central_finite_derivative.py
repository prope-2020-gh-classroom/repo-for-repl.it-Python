from central_finite_derivative import approx_first_derivative, approx_second_derivative
import math

f = math.atan
x = 0.9
h = 1e-6

res_first_d = approx_first_derivative(f,x,h)
res_second_d = approx_second_derivative(f,x,h)
print(res_first_d)
print(res_second_d)