### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 173388ab-207a-42a6-b364-b2c1cb335f6b
# ╠═╡ show_logs = false
begin
	ENV["JULIA_PKG_USE_CLI_GIT"] = true

	import MarkdownLiteral: @mdx

	using Pkg
	using Downloads
	using TOML
	using Distributions
	using PlutoUI
	using Random
	using Base64
	using Plots
	using Test
	
	default(fontfamily="Computer Modern", framestyle=:box) # LaTeX-style plotting

	md"> _Package management._"
end

# ╔═╡ 8f08e006-2d70-4173-a02e-267e8486f5c8
using StanfordAA228V

# ╔═╡ f53306f6-1072-40ae-bf65-1927e5eae088
md"""
# Project 0: Falsification introduction
_A light-weight introduction to falsification._

**Task**: Simply count the number of failures for a 1D Gaussian environment.
- Write a function `num_failures(sys, ψ; m)` that given a system and specification, returns the number of failures over `m` rollouts.
"""

# ╔═╡ 2eec992f-0581-4b7b-a14f-ad119c230ae2
Markdown.MD(HTML("<h3 id='submission'>Submission</h2>"),
Markdown.parse("""
Once completed, you will submit **this file** (`project0.jl`) to Gradescope.

The directory of the notebook is located here:"""),
md"""
- $(@bind directory_trigger OpenDirectory(@__DIR__))
    - ↑ Click to open directory.""",
md"""
If you encounter issues, [please ask us on Ed](https://edstem.org/us/courses/69226/discussion).
""")

# ╔═╡ a3ce3d69-83f0-4ab3-b579-827d793b511e
if directory_trigger
	@info "Opening local directory..."
	sleep(1)
end

# ╔═╡ 17fa8557-9656-4347-9d44-213fd3b635a6
Markdown.parse("""
## System
The system is comprised of an `agent`, environment (`env`), and `sensor`.
""")

# ╔═╡ 22feee3d-4627-4358-9937-3c780b7e8bcb
sys = System(NoAgent(), SimpleGaussian(), IdealSensor())

# ╔═╡ 45f7c3a5-5763-43db-aba8-41ef8db39a53
md"""
## Environment
The environment is a standard normal (Gaussian) distribution $\mathcal{N}(0, 1)$.
"""

# ╔═╡ 9c1daa96-76b2-4a6f-8d0e-f95d26168d2b
ps = Ps(sys.env)

# ╔═╡ ab4c6807-5b4e-4688-b794-159e26a1599b
ψ = LTLSpecification(@formula □(s->s > -2));

# ╔═╡ 370a15eb-df4b-493a-af77-00914b4616ea
Markdown.parse("""
## Specification \$\\psi\$
The specification \$\\psi\$ (written `\\psi<TAB>` in code) indicates what the system should do:

\$\$\\psi(\\tau) = \\square(s > $(ψ.formula.ϕ.c))\$\$

i.e., "the state \$s\$ in the trajectory \$\\tau\$ should _always_ (\$\\square\$) be greater than \$$(ψ.formula.ϕ.c)\$, anything else is a failure."
""")

# ╔═╡ 166bd412-d433-4dc9-b874-7359108c0a8b
Markdown.parse("""
A failure is unlikely given that the probability of failure is:

\$\$P(s > $(ψ.formula.ϕ.c)) \\approx $(round(cdf(ps, ψ.formula.ϕ.c), sigdigits=4))\$\$
""")

# ╔═╡ 00d4d678-a19d-4bba-b8f5-79d7e1466a63
md"""
# Useful interface functions
The following functions are provided by `StanfordAA228V.jl` that you may use.

## `rollout`
**`rollout(sys::System)::Array`** — Run a single rollout of the system `sys` to a depth of `d`.
- `τ` is written as `\tau<TAB>` in code.
- **Note**, the 1D Gaussian only needs to run for a depth of `d=1`.
```julia
function rollout(sys::System; d=1)
    s = rand(Ps(sys.env))
    τ = []
    for t in 1:d
        o, a, s′ = step(sys, s)
        push!(τ, (; s, o, a))
        s = s′
    end
    return τ
end
```

## `isfailure`

**`isfailure(ψ, τ)::Bool`** — Using the specification `ψ`, check if the trajector `τ` led to a failure.
- `ψ` is written as `\psi<TAB>` in code.
```julia
isfailure(ψ::Specification, τ) = !evaluate(ψ, τ)
```
"""

# ╔═╡ ea2f1380-6071-4a20-9b87-daf7d2b7ee33
html_half_space()

# ╔═╡ 86db41bf-c699-426c-a026-971b79dc0e2c
md"""
# 👩‍💻 **Task**: Count the number of failures
Please fill in the following `num_failures` function.
"""

# ╔═╡ ddc8031f-fd06-4189-b00d-70f930998db4
start_code()

# ╔═╡ 4e96d96e-d2c3-4780-8e4d-fbe31503574e
md"""
	num_failures(sys, ψ; m)

A function that takes in a system `sys` and a specification `ψ` and **returns the number of failures**.

- `m` = number of rollouts

**Note**: `ψ` is written as `\psi<TAB>`
"""

# ╔═╡ 798451be-5646-4b5e-b4d7-04d9fc9e6699
function num_failures(sys, ψ; m=1000)
	# TODO: WRITE YOUR CODE HERE
end

# ╔═╡ 651313a4-2766-49dd-8737-475ed80079e2
end_code()

# ╔═╡ c683b6a8-0c0f-4232-b914-a70f3e2b06e8
html_half_space()

# ╔═╡ 873c99d8-ebd8-4ce3-92ca-6975c713fc8b
md"""
## Example usage of `num_failures`
Example usage with `m=1000` number of rollouts.
"""

# ╔═╡ a6e52a4e-6e75-4ae0-9e3a-cc82f9ad6b2b
num_failures(sys, ψ; m=1000)

# ╔═╡ bfc1779c-5ce8-47f5-be65-e4e74f2071cf
hint(md"Can you write the `num_failures` function in one line of code? (Not required)"; title="Feeling adventurous?")

# ╔═╡ eee07ddd-ee38-4a22-8807-3afcf46ad00d
html_quarter_space()

# ╔═╡ 2827a6f3-47b6-4e6f-b6ae-63271715d1f3
Markdown.parse("""
# 📊 Tests
The tests below run your `num_failures` function to see if it works properly.

This will automatically run anytime the `num_failures` function is changed and saved (due to Pluto having dependent cells).
""")

# ╔═╡ 4a91853f-9685-47f3-998a-8e0cfce688f8
Markdown.parse("""
## Running tests
Run two tests, controlling the RNG seed for deterministic outputs.
""")

# ╔═╡ 3edcc3a0-6574-4784-8698-8cec43932935


# ╔═╡ ba6c082b-6e62-42fc-a85c-c8b7efc89b88
# ╠═╡ show_logs = false
begin
	########################################################
	# NOTE: DECODING THIS IS A VIOLATION OF THE HONOR CODE.
	########################################################
	ModuleTA = "UsingThisViolatesTheHonorCode_$(basename(tempname()))"
	try
		eval(Meta.parse("""
		module $ModuleTA
		$(String(base64decode("IyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMKIyBMT09LSU5HIEFUIFRISVMgSVMgQSBWSU9MQVRJT04gT0YgVEhFIEhPTk9SIENPREUKIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMjIyMKClRoaXNNb2R1bGUgPSBzcGxpdChzdHJpbmcoQF9fTU9EVUxFX18pLCAiLiIpW2VuZF0KCiMgTG9hZCBhbGwgY29kZSBhbmQgcGFja2FnZXMgZnJvbSBwYXJlbnQgbW9kdWxlClBhcmVudCA9IHBhcmVudG1vZHVsZShAX19NT0RVTEVfXykKCm1vZHVsZXMobTo6TW9kdWxlKSA9IGNjYWxsKDpqbF9tb2R1bGVfdXNpbmdzLCBBbnksIChBbnksKSwgbSkKCiMgTG9hZCBmdW5jdGlvbnMgYW5kIHZhcmlhYmxlcwpmb3IgbmFtZSBpbiBuYW1lcyhQYXJlbnQsIGltcG9ydGVkPXRydWUpCglpZiBuYW1lICE9IFN5bWJvbChUaGlzTW9kdWxlKSAmJiAhb2NjdXJzaW4oIiMiLCBzdHJpbmcobmFtZSkpICYmICFvY2N1cnNpbigiVXNpbmdUaGlzVmlvbGF0ZXNUaGVIb25vckNvZGUiLCBzdHJpbmcobmFtZSkpCgkJQGV2YWwgY29uc3QgJChuYW1lKSA9ICQoUGFyZW50KS4kKG5hbWUpCgllbmQKZW5kCgpleGNsdWRlcyA9IFsiUGx1dG9SdW5uZXIiLCAiSW50ZXJhY3RpdmVVdGlscyIsICJNYXJrZG93biIsICJDb3JlIiwgIkJhc2UiLCAiQmFzZS5NYWluSW5jbHVkZSJdCgojIExvYWQgcGFja2FnZXMKZm9yIG1vZCBpbiBtb2R1bGVzKFBhcmVudCkKCXN0cmluZyhtb2QpIGluIGV4Y2x1ZGVzICYmIGNvbnRpbnVlCgl0cnkKCQlAZXZhbCB1c2luZyAkKFN5bWJvbChtb2QpKQoJY2F0Y2ggZXJyCgkJaWYgZXJyIGlzYSBBcmd1bWVudEVycm9yCgkJCXRyeQoJCQkJQGV2YWwgdXNpbmcgU3RhbmZvcmRBQTIyOFYuJChTeW1ib2wobW9kKSkKCQkJY2F0Y2ggZXJyMgoJCQkJQHdhcm4gZXJyMgoJCQllbmQKCQllbHNlCgkJCUB3YXJuIGVycgoJCWVuZAoJZW5kCmVuZAoKdXNpbmcgU3RhbmZvcmRBQTIyOFYKCmZ1bmN0aW9uIG51bV9mYWlsdXJlc190ZXN0KHN5cywgz4g7IG09MTAwMCkKICAgIHJldHVybiBzdW0oaXNmYWlsdXJlLijPiCwgcm9sbG91dChzeXMpIGZvciBpIGluIDE6bSkpCmVuZAo=")))
		end"""))
		global UsingThisViolatesTheHonorCode = getfield(@__MODULE__, Symbol(ModuleTA))
	catch err
		@warn err
	end
	𝑡𝑟𝑖𝑔𝑔𝑒𝑟 = true
	global SEED = sum(Int.(collect("AA228V"))) # Cheeky seed value :)
	md"""
	# Backend
	_Helper functions and project management._
	"""
end

# ╔═╡ 83884eb4-6718-455c-b731-342471325326
function run_project0_test(num_failures::Function; m=1000, seed=SEED)
	Random.seed!(seed) # For determinism
	return num_failures(sys, ψ; m)
end

# ╔═╡ b6f15d9c-33b8-40e3-be57-d91eda1c9753
begin
	test1_m = 1000
	test1_output = run_project0_test(num_failures; m=test1_m, seed=SEED)
end

# ╔═╡ 2ff6bb9c-5282-4ba1-b62e-a9fd0fe1969c
Markdown.parse("""
### Test 1: \$m = $test1_m\$
""")

# ╔═╡ 522bb285-bc06-4c92-82ee-c0d0f68b184c
if isa(test1_output, Number)
	Markdown.parse("""
	The frequentist failure probability estimate for test 1 would be:
	
	\$\$\\hat{p}_{\\rm failure}^{({\\rm test}_1)} = \\frac{$(test1_output)}{$test1_m} =  $(round(test1_output/test1_m, sigdigits=5))\$\$
	""")
else
	md"*Update `num_failures` to get an estimated failure probability for test 1.*"
end

# ╔═╡ 3314f402-10cc-434c-acbc-d38e59e4b846
begin
	test2_m = 5000
	test2_output = run_project0_test(num_failures; m=test2_m, seed=SEED)
end

# ╔═╡ 089581ec-8aff-4c56-9a65-26d394d5eec3
Markdown.parse("""
### Test 2: \$m = $test2_m\$
""")

# ╔═╡ d72be566-6ad7-4817-8590-a504a699a4da
if isa(test2_output, Number)
	Markdown.parse("""
	The frequentist failure probability estimate for test 2 would be:
	
	\$\$\\hat{p}_{\\rm failure}^{({\\rm test}_2)} = \\frac{$(test2_output)}{$test2_m} =  $(round(test2_output/test2_m, sigdigits=5))\$\$
	""")
else
	md"*Update `num_failures` to get an estimated failure probability for test 2.*"
end

# ╔═╡ 956603a8-5a5f-41b8-ba52-b0a678433342
𝑓 = UsingThisViolatesTheHonorCode.num_failures_test; 𝑡𝑟𝑖𝑔𝑔𝑒𝑟;

# ╔═╡ 6302729f-b34a-4a18-921b-d194fe834208
begin
	# ⚠️ Note: PLEASE DO NOT MODIFY. Why are you in here anyhow :)?

	if isnothing(test1_output) && isnothing(test1_output)
		test1_passed = test2_passed = false
		almost()
	else
		test1_passed = test1_output == run_project0_test(𝑓; m=test1_m, seed=SEED)
		test2_passed = test2_output == run_project0_test(𝑓; m=test2_m, seed=SEED)
		if test1_passed && test2_passed
			correct(md"""
    All tests have passed, you're done with Project 0!

    Please submit `project0.jl` (this file) to Gradescope (see the [Submission](#submission) section).
    """)
		else
			keep_working()
		end
	end
end

# ╔═╡ cee165f0-049f-4ea3-8f19-04e66947a397
HTML("""
<h3>$(test1_passed && test2_passed ? "✔️" : "❌") Check tests</h3>
<p>If the following test indicator is <span style='color:#759466'><b>green</b></span>, you can submit <code>project0.jl</code> (this file) to Gradescope.</p>
""")

# ╔═╡ 8461e634-8830-477c-99c8-894dea63d5ce
begin
	if !ismissing(directory_trigger) && directory_trigger
		try
			if Sys.iswindows()
				run(`explorer $(abspath(@__DIR__))`)
			elseif Sys.isapple()
				run(`open $(abspath(@__DIR__))`)
			elseif Sys.islinux()
				run(`xdg-open $(abspath(@__DIR__))`)
			end
		catch end
	end
end; md"> _Helper for opening local directories._"

# ╔═╡ ef084fea-bf4d-48d9-9c84-8cc1dd98f2d7
Markdown.MD(TableOfContents(), md"> _Table of contents._")

# ╔═╡ a6931d1e-08ad-4592-a54c-fd76cdc51294
Markdown.MD(md"$(@bind dark_mode DarkModeIndicator())", md"> _Dark mode indicator._")

# ╔═╡ 28504aea-ee04-4b70-a067-a450d8d6374f
try
	if !validate_version(StanfordAA228V)
		Markdown.MD(
			almost(md"""
			Your `StanfordAA228V` package is out-of-date. Please update it via the instructions below.

			**Then restart the notebook.**

			_(This warning may persist after restart, wait until the notebook finishes loading entirely)_"""),				md"""$(LocalResource(joinpath(@__DIR__, "..", "media", dark_mode ? "update-package-dark.gif" : "update-package.gif")))"""
		)
	end
catch end

# ╔═╡ 0cdadb29-9fcd-4a70-9937-c24f07ce4657
begin
	if dark_mode
		plot(
			bg="transparent",
			background_color_inside="black",
			bglegend="black",
			fg="white",
			gridalpha=0.5,
		)
	else
		plot()
	end

	# Create a range of x values
	_X = range(-4, 4, length=1000)
	_Y = pdf.(ps, _X)
	
	# Plot the Gaussian density
	plot!(_X, _Y,
	     xlim=(-4, 4),
	     ylim=(-0.001, 0.41),
	     linecolor=dark_mode ? "white" : "black",
		 fillcolor=dark_mode ? "darkgray" : "lightgray",
		 fill=true,
	     xlabel="state \$s\$",
	     ylabel="density \$p(s)\$",
	     size=(600, 300),
	     label=false)
	
	# Identify the indices where x <= -2
	idx = _X .<= ψ.formula.ϕ.c
	
	# Extract the x and y values for the region to fill
	x_fill = _X[idx]
	y_fill = _Y[idx]
	
	# Create the coordinates for the filled polygon
	# Start with the x and y values where x <= -2
	# Then add the same x values in reverse with y = 0 to close the polygon
	polygon_x = vcat(x_fill, reverse(x_fill))
	polygon_y = vcat(y_fill, zeros(length(y_fill)))
	
	# Add the filled area to the plot
	plot!(polygon_x, polygon_y,
	      fill=true,
	      fillcolor="crimson",
	      linecolor="transparent", # No border for the filled area
		  alpha=0.5,
	      label=false)

	# Draw failure threshold
	vline!([ψ.formula.ϕ.c], color="crimson", label="Failure threshold")
end

# ╔═╡ c9c45286-58a4-40e6-b2a4-d828e627c6ec
Markdown.MD(html"""
<style>
	h3 {
		border-bottom: 1px dotted var(--rule-color);
	}

    .container {
      display: flex;
      align-items: center;
      width: 100%;
      margin: 1px 0;
    }

    .line {
      flex: 1;
      height: 2px;
      background-color: #B83A4B;
    }

    .text {
      margin: 0 5px;
      white-space: nowrap; /* Prevents text from wrapping */
    }
</style>""", md"> _Notebook styling._")

# ╔═╡ Cell order:
# ╟─28504aea-ee04-4b70-a067-a450d8d6374f
# ╠═8f08e006-2d70-4173-a02e-267e8486f5c8
# ╟─f53306f6-1072-40ae-bf65-1927e5eae088
# ╟─2eec992f-0581-4b7b-a14f-ad119c230ae2
# ╟─a3ce3d69-83f0-4ab3-b579-827d793b511e
# ╟─17fa8557-9656-4347-9d44-213fd3b635a6
# ╠═22feee3d-4627-4358-9937-3c780b7e8bcb
# ╟─45f7c3a5-5763-43db-aba8-41ef8db39a53
# ╠═9c1daa96-76b2-4a6f-8d0e-f95d26168d2b
# ╟─370a15eb-df4b-493a-af77-00914b4616ea
# ╠═ab4c6807-5b4e-4688-b794-159e26a1599b
# ╟─0cdadb29-9fcd-4a70-9937-c24f07ce4657
# ╟─166bd412-d433-4dc9-b874-7359108c0a8b
# ╟─00d4d678-a19d-4bba-b8f5-79d7e1466a63
# ╟─ea2f1380-6071-4a20-9b87-daf7d2b7ee33
# ╟─86db41bf-c699-426c-a026-971b79dc0e2c
# ╟─ddc8031f-fd06-4189-b00d-70f930998db4
# ╟─4e96d96e-d2c3-4780-8e4d-fbe31503574e
# ╠═798451be-5646-4b5e-b4d7-04d9fc9e6699
# ╟─651313a4-2766-49dd-8737-475ed80079e2
# ╟─c683b6a8-0c0f-4232-b914-a70f3e2b06e8
# ╟─873c99d8-ebd8-4ce3-92ca-6975c713fc8b
# ╠═a6e52a4e-6e75-4ae0-9e3a-cc82f9ad6b2b
# ╟─bfc1779c-5ce8-47f5-be65-e4e74f2071cf
# ╟─eee07ddd-ee38-4a22-8807-3afcf46ad00d
# ╟─2827a6f3-47b6-4e6f-b6ae-63271715d1f3
# ╠═83884eb4-6718-455c-b731-342471325326
# ╟─4a91853f-9685-47f3-998a-8e0cfce688f8
# ╟─2ff6bb9c-5282-4ba1-b62e-a9fd0fe1969c
# ╠═b6f15d9c-33b8-40e3-be57-d91eda1c9753
# ╟─522bb285-bc06-4c92-82ee-c0d0f68b184c
# ╟─089581ec-8aff-4c56-9a65-26d394d5eec3
# ╠═3314f402-10cc-434c-acbc-d38e59e4b846
# ╟─d72be566-6ad7-4817-8590-a504a699a4da
# ╟─3edcc3a0-6574-4784-8698-8cec43932935
# ╟─cee165f0-049f-4ea3-8f19-04e66947a397
# ╟─6302729f-b34a-4a18-921b-d194fe834208
# ╟─956603a8-5a5f-41b8-ba52-b0a678433342
# ╟─ba6c082b-6e62-42fc-a85c-c8b7efc89b88
# ╟─173388ab-207a-42a6-b364-b2c1cb335f6b
# ╟─8461e634-8830-477c-99c8-894dea63d5ce
# ╟─ef084fea-bf4d-48d9-9c84-8cc1dd98f2d7
# ╟─a6931d1e-08ad-4592-a54c-fd76cdc51294
# ╟─c9c45286-58a4-40e6-b2a4-d828e627c6ec
