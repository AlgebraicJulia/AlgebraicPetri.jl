### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 86ffd357-1510-4d05-8a38-b59b42b79b39
begin
using Pkg
Pkg.activate(".")
end

# ╔═╡ 3779b846-e5ec-4239-a1d4-af2f8c2f10eb
begin
	include("EnzymeReactions.jl")
	import .EnzymeReactions: ob, ode,
							 inactivate, bindunbind, degrade,
							 enzX, enzXY, enzXsubY,
							 enz, enz_enz, enz_sub,
							 enzyme_uwd
	using PlutoUI

	using AlgebraicPetri
	using Catlab.WiringDiagrams
	using Catlab.CategoricalAlgebra
	using Catlab.Graphics
	using LabelledArrays

	using DifferentialEquations
	using Plots
	using Interact

	display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));
	nothing
end

# ╔═╡ 93df89f0-8429-4fcc-bd01-6982417f5134
begin
	# Initial Concentrations
	K = :K=>33000;
	S = :S=>33000;
	L = :L=>33000;
	Kinact = :Kinact=>0;
	Sinact = :Sinact=>0;
	Linact = :Linact=>0;
	E = :E=>700000;
	G = :G=>1300000;
	
	# Parameter Rates (units of pM and min)
	rxns = Dict(
	  :K => [inactivate(K, 7.494e-10)
	         bindunbind(K, K, 7.814e-4, 3.867e-3)
	         degrade(K, K, 2.265e-1)
	         bindunbind(K, Kinact, 7.814e-4, 3.867e-3)
	         degrade(K, Kinact, 2.265e-1)],
	  :S => [inactivate(S, 1.906e-2)
	         bindunbind(S, S, 3.534e-7, 1.688e2)
	         degrade(S, S, 1.433e1)
	         bindunbind(S, Sinact, 3.534e-7, 1.688e2)
	         degrade(S, Sinact, 1.433e1)],
	  :L => [inactivate(L, 7.810e-3)
	         bindunbind(L, L, 1.000e-11, 7.440e3)
	         degrade(L, L, 2.670e2)
	         bindunbind(L, Linact, 1.000e-11, 7.440e3)
	         degrade(L, Linact, 2.670e2)],
	  :KE => [bindunbind(K, E, 9.668e-6,1.000e-2)
	          degrade(K, E, 1.728e0)],
	  :KG => [bindunbind(K, G, 2.764e-6, 8.780e-1)
	          degrade(K, G, 1.502e0)],
	  :SE => [bindunbind(S, E, 4.197e-7, 1.06e-3)
	          degrade(S, E, 1.384e4)],
	  :SG => [bindunbind(S, G, 5.152e-8, 3.894e-3)
	          degrade(S, G, 8.755e-1)],
	  :LE => [bindunbind(L, E, 1.977e-8, 1.000e-2)
	          degrade(L, E, 1.066e2)],
	  :LG => [bindunbind(L, G, 3.394e-8, 2.365e1)
	          degrade(L, G, 4.352e0)],
	  :KS => [bindunbind(K, S, 8.822e-4, 4.114e5)
	          degrade(K, S, 9.000e-10)
	          bindunbind(K, Sinact, 8.822e-4, 4.114e5)
	          degrade(K, Sinact, 9.000e-10)],
	  :KL => [bindunbind(K, L, 1.756e-4, 3.729e4)
	          degrade(K, L, 6.505e6)
	          bindunbind(K, Linact, 1.756e-4, 3.729e4)
	          degrade(K, Linact, 6.505e6)],
	  :SK => [bindunbind(S, K, 3.679e-4, 1.562e3)
	          degrade(S, K, 4.410e2)
	          bindunbind(S, Kinact, 3.679e-4, 1.562e3)
	          degrade(S, Kinact, 4.410e2)],
	  :SL => [bindunbind(S, L, 1.000e-3, 5.000e2)
	          degrade(S, L, 1.000e-7)
	          bindunbind(S, Linact, 1.000e-3, 5.000e2)
	          degrade(S, Linact, 1.000e-7)],
	  :LK => [bindunbind(L, K, 1.000e-3, 4.118e3)
	          degrade(L, K, 3.234e1)
	          bindunbind(L, Kinact, 1.000e-3, 4.118e3)
	          degrade(L, Kinact, 3.234e1)],
	  :LS => [bindunbind(L, S, 1.056e-12, 5.000e2)
	          degrade(L, S, 5.000e-1)
	          bindunbind(L, Sinact, 1.056e-12, 5.000e2)
	          degrade(L, Sinact, 5.000e-1)]
	);
	
	# define labels to reaction network mappings
	functor(x) = oapply(x, Dict(
	  :catK=>enz(rxns, K),
	  :catS=>enz(rxns, S),
	  :catL=>enz(rxns, L),
	  :catKcatS=>enz_enz(rxns, K,S),
	  :catKcatL=>enz_enz(rxns, K,L),
	  :catScatK=>enz_enz(rxns, S,K),
	  :catScatL=>enz_enz(rxns, S,L),
	  :catLcatK=>enz_enz(rxns, L,K),
	  :catLcatS=>enz_enz(rxns, L,S),
	  :catKsubE=>enz_sub(rxns, K,E),
	  :catSsubE=>enz_sub(rxns, S,E),
	  :catLsubE=>enz_sub(rxns, L,E),
	  :catKsubG=>enz_sub(rxns, K,G),
	  :catSsubG=>enz_sub(rxns, S,G),
	  :catLsubG=>enz_sub(rxns, L,G)));
	nothing
end

# ╔═╡ e6589d31-dce7-42c3-b494-db03fe561ae9
	uwd = enzyme_uwd([:K, :S], [:E]);

# ╔═╡ 7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
begin
	display_uwd(uwd)
end

# ╔═╡ cf9e03db-42b7-41f6-80ce-4b12ddb93211
begin
	model = uwd |> functor |> apex;
	r = join(["'$k': $(rates(model)[k])" for k in keys(rates(model))], ", ");
	#== Basic Slider Template
<li><input
	class='slider'
	type='range' 
	min='-10'
	step='0.01' 
	max='10'
	value='0'
	oninput='this.nextElementSibling.value=(10**this.value).toExponential(2)'
	></input>
<output>1</output></li>
==#

form_vals = HTML("
<form>
	<div style='height:500px;overflow:scroll'>
	<table>
	</table>
	</div>
</form>
<script>
//`currentScript` is the current script tag - we use it to select elements//
const form = currentScript.parentElement.querySelector('form')
const list = form.querySelector('table')
var rates = {$r};
console.log(rates)
for ( var r in rates ){
	console.log(r);
	var item = document.createElement('tr');
	var label = document.createElement('th');
	label.innerText = r
	var slider_box = document.createElement('td');
	var slider_val_box = document.createElement('td');
	var slider = document.createElement('input');
	var slider_val = document.createElement('input');
	slider_val.setAttribute('type', 'text');
	slider_val.setAttribute('class', 'sliderval');
	slider_val.value = rates[r].toExponential(2);
	slider.setAttribute('type', 'range');
	slider.setAttribute('class', 'slider');
	slider.setAttribute('min', '-9.99');
	slider.setAttribute('max', '10');
	slider.setAttribute('value', '0.0');
	slider.setAttribute('step', '0.01');
	slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${rates[r]}).toExponential(2)`);
	slider_box.appendChild(slider)
	slider_val_box.appendChild(slider_val)
	item.appendChild(label)
	item.appendChild(slider_box)
	item.appendChild(slider_val_box)
	list.appendChild(item)
}
var x = form.getElementsByClassName('sliderval');
function onsubmit(){
	// We send the value back to Julia //
	form.value = Array.from(x, (v,_)=>{return v.value})
	form.dispatchEvent(new CustomEvent('input'))
	console.log(form.value)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.value = 'Update Plot';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()
</script>
");
nothing
end

# ╔═╡ ba87cd7e-e9c7-4a20-99be-eee794f968a1
@bind c form_vals

# ╔═╡ 066b7505-e21b-467e-86c1-cea1ff80246e
begin
	cur_rate = Dict(tnames(model)[i]=>parse(Float64, c[i]) for i in 1:length(c))
	cur_conc = concentrations(model)
  	vf = vectorfield(model);
	nothing
end

# ╔═╡ 1ba7bbe5-7a85-454e-a9cf-deaf5f00d6ad
sol = solve(ODEProblem(vf, cur_conc, (0.0,120.0),cur_rate));

# ╔═╡ a141cd27-6ea0-4f73-80b5-72d8e5770ed4
begin
	tsteps = range(0.0,120.0, length=5000)
	labels = [:K, :Kdeg, :KK]
	plot([[sol(t)[l] for t in tsteps] for l in labels], labels=hcat(String.(labels)...), linewidth=3)
end

# ╔═╡ Cell order:
# ╟─86ffd357-1510-4d05-8a38-b59b42b79b39
# ╟─3779b846-e5ec-4239-a1d4-af2f8c2f10eb
# ╟─93df89f0-8429-4fcc-bd01-6982417f5134
# ╠═e6589d31-dce7-42c3-b494-db03fe561ae9
# ╟─7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
# ╟─cf9e03db-42b7-41f6-80ce-4b12ddb93211
# ╠═ba87cd7e-e9c7-4a20-99be-eee794f968a1
# ╟─066b7505-e21b-467e-86c1-cea1ff80246e
# ╟─1ba7bbe5-7a85-454e-a9cf-deaf5f00d6ad
# ╠═a141cd27-6ea0-4f73-80b5-72d8e5770ed4
