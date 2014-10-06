function x=minabs(a,b)
x=0.5*(sign(abs(b)-abs(a))*(a-b)+a+b);