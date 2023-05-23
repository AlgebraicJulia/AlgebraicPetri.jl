module AlgebraicPetriOrdinaryDiffEqExt

using AlgebraicPetri
import OrdinaryDiffEq: ODEProblem

ODEProblem(m::AbstractPetriNet, u0, tspan, β) = ODEProblem(vectorfield(m), u0, tspan, β)
ODEProblem(m::AbstractPetriNet, tspan) = ODEProblem(vectorfield(m), concentrations(m), tspan, rates(m))

end
