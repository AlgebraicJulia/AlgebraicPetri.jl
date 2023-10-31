module AlgebraicPetriOrdinaryDiffEqExt

using AlgebraicPetri
import OrdinaryDiffEq: ODEProblem

""" ODEProblem(m::AbstractPetriNet, u0, tspan, β)

Convert a general PetriNet to an ODEProblem
This takes the initial conditions and rates as parameters and will ignore
any initial conditions and rates embedded in the Petri Net if they exist.
"""
ODEProblem(m::AbstractPetriNet, u0, tspan, β) = ODEProblem(vectorfield(m), u0, tspan, β)

""" ODEProblem(m::AbstractPetriNet, tspan)

Convert a general PetriNet to an ODEProblem
This assumes the concentrations and weights are embedded in the Petri Net.
"""
ODEProblem(m::AbstractPetriNet, tspan) = ODEProblem(vectorfield(m), concentrations(m), tspan, rates(m))

end
