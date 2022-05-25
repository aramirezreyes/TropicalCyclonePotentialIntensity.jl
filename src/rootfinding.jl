"""
    approximate_t_newton_raphson(tguess,penv,target_entropy)
Given a parcel's total entropy s_parcel and a pressure level P, use the Newton-Raphson method to find T⋆ such that s(T⋆,P) = s_parcel. It assumes reversible moist process.
Total entropy of the parcel before condensation: 
s = (cpd + rt*cl)ln(T) - Rd*ln(pd) + Lv*r/T - r*Rvln(H) : rt is total water (vapor plus liquid) 
"""
function find_root_newton_raphson(func, func_derivative; target_value = 0.0,  initial_guess = 0.0, atol = 0.001*unit(initial_guess), max_iter = 500)
    niter = 0
    initial_err = atol + 5unit(initial_guess)
    approximation = initial_guess
    err = initial_err
    step_size = 1.0
    while (abs(err) > 0.001 * unit(initial_guess))
        niter += 1
        step_size = niter < 3 ? 0.3 : 1.0
        approximation = approximation - step_size*(func(approximation) - target_value)/func_derivative(approximation)
        err = step_size*(func(approximation) - target_value)/func_derivative(approximation)
        if (niter > max_iter ) 
            error("Function didn't converge after $max_iter iterations")
        end
    end
    return approximation
end
