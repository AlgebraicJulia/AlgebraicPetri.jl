module solvers

export ODEProblem

using DiffEqBase
using OrdinaryDiffEq
using SteadyStateDiffEq
using StochasticDiffEq
using JumpProcesses
using SparseArrays
using AlgebraicPetri

import OrdinaryDiffEq: ODEProblem
import SteadyStateDiffEq: SteadyStateProblem
import StochasticDiffEq: SDEProblem
import JumpProcesses: JumpProblem

# Helper function to handle both constant and time-dependent rates
valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end

""" ODEProblem(pn::AbstractPetriNet, u0, tspan, β)

Generate an OrdinaryDiffEq ODEProblem from a Petri net
*Note*: This currently exists in another ext module: AlgebraicPetriOrdinaryDiffEqExt, 
        so it is being omitted for now
"""
# ODEProblem(pn::AbstractPetriNet, u0, tspan, β) = ODEProblem(vectorfield(pn), u0, tspan, β)

# Continuous state-dependent callback function to reset states to zero
function statecb(s)
    # Condition: when the state u[s] reaches zero
    cond = (u, t, integrator) -> u[s]
    # Affect: set u[s] to zero
    aff = (integrator) -> integrator.u[s] = 0.0
    return ContinuousCallback(cond, aff)
end

""" SteadyStateProblem(pn::AbstractPetriNet, u0, tspan, β)

Generate an SteadyStateDiffEq SteadyStateProblem from a Petri net
"""
SteadyStateProblem(pn::AbstractPetriNet, u0, tspan, β) = SteadyStateProblem(ODEProblem(vectorfield(pn), u0, tspan, β))


""" SDEProblem(pn::AbstractPetriNet, u0, tspan, β)

Generate an StochasticDiffEq SDEProblem and an appropiate CallbackSet
"""

function SDEProblem(pn::AbstractPetriNet, u0, tspan, β)
    tm = TransitionMatrices(pn)
    input = tm.input
    output = tm.output

    num_transitions, num_states = size(input)
    nu = spzeros(Float64, num_states, num_transitions)  # Noise matrix (state x transition)

    # Set up the matrix `nu` for the noise term
    for tr in 1:num_transitions  
        for st in 1:num_states  
            nu[st, tr] -= input[tr, st]
            nu[st, tr] += output[tr, st]
        end
    end

    # Noise function
    noise(du, u, p, t) = begin

        rates = zeros(valtype(du), num_transitions)
        u_m = [u[sname(pn, i)] for i in 1:num_states]
        p_m = [p[tname(pn, i)] for i in 1:num_transitions]
        for tr in 1:num_transitions
            rates[tr] = valueat(p_m[tr], u, t) * prod(u_m[st]^input[tr, st] for st in 1:num_states)
        end

        # Update du
        for tr in 1:num_transitions
            rate_sqrt = sqrt(abs(rates[tr]))  # sqrt of the absolute transition rate
            for st in 1:num_states
                du[st, tr] = rate_sqrt * (output[tr, st] - input[tr, st])
            end
        end

        return du
    end

    # Create a CallbackSet for all states
    callback_set = CallbackSet([statecb(s) for s in 1:num_states]...)  # Apply statecb for each state `s`

    # Return the SDEProblem with the vector field, noise, and callbacks
    return SDEProblem(vectorfield(pn), noise, u0, tspan, β, noise_rate_prototype=nu),
    callback_set
    
end

# Calculate the rate for a single transition
function transition_rate(pn::AbstractPetriNet, tm::TransitionMatrices, u, p, t, i)
    # Extract the state (u_m) and parameter (p_m) for the Petri net
    u_m = [u[sname(pn, j)] for j in 1:ns(pn)]  # Current states
    p_m = [p[tname(pn, j)] for j in 1:nt(pn)]  # All transitions

    # Compute the rate for transition `i`
    rate = valueat(p_m[i], u, t) * prod(u_m[j]^tm.input[i, j] for j in 1:ns(pn))
    
    return rate
end

# Wrapper function for jump transition rates
function jumpTransitionRate(pn::AbstractPetriNet, tm::TransitionMatrices, tr)
    return (u, p, t) -> transition_rate(pn, tm, u, p, t, tr)
end

# Function to update the state after a jump occurs
function jumpTransitionFunction(input, output, tr)
    return (integrator) -> begin
        for st in 1:length(integrator.u)
            integrator.u[st] -= input[tr, st]
            integrator.u[st] += output[tr, st]
        end
    end
end


""" JumpProblem(pn::AbstractPetriNet, u0, tspan, β)

Generate an JumpProcesses JumpProblem from a Petri net
"""
function JumpProblem(pn::AbstractPetriNet, u0, tspan, p)
    tm = TransitionMatrices(pn)
    num_transitions = nt(pn)
    input = tm.input
    output = tm.output
    prob = DiscreteProblem(u0, tspan, p)
    return JumpProblem(prob, Direct(), [ConstantRateJump(jumpTransitionRate(pn, tm, tr), 
                                        jumpTransitionFunction(input, output, tr)) for tr in 1:num_transitions]...)

end

end