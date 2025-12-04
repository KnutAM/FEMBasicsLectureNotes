### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 6b0414f5-13f1-4e68-944d-d29b18affaa3
begin
	import CairoMakie as Plt
	import GeometryBasics as GB
	using PlutoUI
	using PlutoTeachingTools
	using Ferrite
	using FerriteAssembly
	using MechanicalMaterialModels
	using MaterialModelsBase
	using LaTeXStrings
	using Printf
	using Format
end

# ╔═╡ 13c0d238-e1cc-473c-a41f-dfbb05622b0f
begin
	# Some plotting defaults (works best at the beginning)
	fontsize_theme = Plt.merge(Plt.theme_latexfonts(), Plt.Theme(fontsize = 18))
	elemsize_theme = Plt.Theme(;linewidth = 3, markersize = 12)
	Plt.set_theme!(Plt.merge(fontsize_theme, elemsize_theme))
	# Display navigation
	TableOfContents()
end

# ╔═╡ b92cd124-af25-11f0-b55d-45d572d7ae78
md"""
# Welcome to FEM Basics (VSM167)
This material is provided under the MIT License: [license file](https://raw.githubusercontent.com/KnutAM/FEMBasicsLectureNotes/refs/heads/main/LICENSE). If you find any errors, please [open an issue on github](https://github.com/KnutAM/FEMBasicsLectureNotes/issues/new), or [send me an email](mailto:knut.andreas.meyer@chalmers.se).
"""

# ╔═╡ 49a9b55e-16ad-4d80-ac5d-cbf18cf717d7
md"""
## L1: Reference element

### Shape functions
The Finite Element Method (FEM) is all about approximating the solution to Partial Differential Equations (PDEs) with a linear combination of shape functions, ``\hat{N}_i(\xi)``, i.e.
```math
       g(\xi) \approx f(\xi; \underline{a}) = \sum_{i = 1}^{\hat{N}_\mathrm{s}} N_i(\xi) a_i
```
Here, ``\hat{N}_i(\xi)`` is the $i$th shape function and ``\underline{a} = [a_1, a_2, \cdots, a_{N_\mathrm{s}}]^\mathrm{T}`` is the coefficient (Degree of Freedom, DoF) vector. For example, we can use two linear shape functions,
```math
\hat{N}_1(\xi) = \frac{1 - \xi}{2},\quad \hat{N}_2(\xi) = \frac{1 + \xi}{2}
```
"""

# ╔═╡ 9bda1921-1001-4bcb-9d53-101fe7010b72
begin
	fig11 = Plt.Figure(size=(500,250))
	ax11 = Plt.Axis(fig11[1,1]; xlabel = L"\xi", ylabel = L"\hat{N}_i(\xi)")
	Plt.lines!(ax11, [-1, 1], [1, 0]; label = L"\hat{N}_1(\xi)")
	Plt.lines!(ax11, [-1, 1], [0, 1]; label = L"\hat{N}_2(\xi)")
	Plt.axislegend(ax11; position=:lc, )
	fig11
end

# ╔═╡ 42924ddb-99ac-4d9d-8092-8b43b8f6b9f4
md"""
We can equally well approximate using nonlinear shape functions, e.g. 
```math
\hat{N}_1(\xi) = \frac{\xi(\xi - 1)}{2}, \quad 
\hat{N}_2(\xi) = \frac{\xi(\xi + 1)}{2}, \quad 
\hat{N}_3(\xi) = 1 - \xi^2
```
"""

# ╔═╡ 6a1731e3-c2bd-41af-8f7b-cc515d2e7912
md"""
Notice that in all cases, we have that each DoF corresponds to the approximated value at a specific point. For the case with two linear shape functions, we have that,
```math
\begin{align}
f(-1; \underline{a}) &= a_1 \\
f(1; \underline{a}) &= a_2
\end{align}
```
These points, (-1) and (+1), are called nodes. For the 3 quadratic shape functions, we have three nodes, (-1), (+1), and (0), corresponding to ``a_1``, ``a_2``, and ``a_3``, such that we have,
```math
\begin{align}
f(-1; \underline{a}) &= a_1 \\
f(1; \underline{a}) &= a_2 \\
f(0; \underline{a}) &= a_3
\end{align}
```
"""

# ╔═╡ 609a415e-0734-4061-997a-4cfb8b71a12f
md"""
### Approximating known functions
For this first step, we will not consider the solution to a PDE, but simply use a linear combination of shape functions to approximate known functions. Specifically, we want to minimize the error,
```math
    E(\underline{a}) = \int_{-1}^{1} \frac{1}{2}\left[f(\xi; \underline{a}) - g(\xi)\right]^2\ \mathrm{d}\xi
```
for different functions $g(\xi)$. In FEM, we choose the shape functions first, and then look for the coefficients that solves our problem. Therefore, we look for the stationary point,
```math
    0 = \frac{\partial E}{\partial a_i} = \underbrace{\left[\int_{-1}^{1} \hat{N}_i(\xi) \hat{N}_j(\xi)\ \mathrm{d}\xi\right]}_{K_{ij}}\ a_j - 
    \underbrace{\int_{-1}^{1} \hat{N}_i(\xi) g(\xi)\ \mathrm{d}\xi}_{f_i}
```
Resulting in the linear equation system
```math
\begin{align}
    \underline{\underline{K}}\ \underline{a} = \underline{f}, \quad \text{or in index notation}^+, \quad K_{ij} a_j = f_i \\
    K_{ij} = \int_{-1}^{1} \hat{N}_i(\xi) \hat{N}_j(\xi)\ \mathrm{d}\xi, \quad 
    f_i = \int_{-1}^{1} \hat{N}_i(\xi) g(\xi)\ \mathrm{d}\xi
\end{align}
```
"""

# ╔═╡ 216d9ced-4729-40a0-9a45-702b85dd9c3d
md"""
!!! note "⁺Index notation"
	The index form is a shorthand to express the indices of the vectors and matrices, such that
	```math
	\underline{f} = \begin{bmatrix} f_1 \\ f_2 \\ \vdots \\ f_N \end{bmatrix},
	\quad
	\underline{\underline{K}} = \begin{bmatrix}
	K_{11} & K_{12} & \cdots & K_{1N} \\
	K_{21} & K_{22} & \cdots & K_{2N} \\
	\vdots & \vdots & \ddots & \vdots \\
	K_{N1} & K_{N2} & \cdots & K_{NN}
	\end{bmatrix}
	```
	And repeated indices in the same term imply summation, i.e. 
	```math
	K_{i\textcolor{red}{j}}a_\textcolor{red}{j} = f_i \Rightarrow K_{i\textcolor{red}{1}}a_\textcolor{red}{1} + K_{i\textcolor{red}{2}}a_\textcolor{red}{2} + \cdots K_{i\textcolor{red}{N}} a_\textcolor{red}{N} = f_i
	```
	where ``\textcolor{red}{j}`` is the repeated index. In this case, that corresponds to the matrix-vector multiplication, ``\underline{\underline{K}}\ \underline{a} = \underline{f}``.

"""

# ╔═╡ 88e934f1-2a7f-478c-8876-f7edbfa1088b
begin
	function setup_numint_plot()
		fig = Plt.Figure()
		title = Plt.Observable(L"title")
		ax = Plt.Axis(fig[1,1]; xlabel = "ξ", ylabel = "h(ξ)", title)
		Plt.xlims!(ax, (-1.2, 1.2))
		ξv = range(-1, 1, 60)
		color = (Plt.Makie.wong_colors()[1], 0.4)
		barplot = Plt.barplot!(ax, [NaN], [NaN]; width = [1.0], gap = 0.0, color, strokewidth=2)
		lineplot = Plt.lines!(ax, ξv, ξv; color = Plt.Makie.wong_colors()[2])
		scatterplot = Plt.scatter!(ax, [NaN], [NaN])
		#textplot = Plt.text!(ax, 0.0, 2.0; text = "N/A", align = (:center, :bottom))
		return fig, (;ξv, barplot, lineplot, scatterplot, title)
	end
	fig_numint, data_numint = setup_numint_plot();
	numint_fun_selector = @bind numint_fun Select([
		(x -> (x - 0.5)) => "(ξ - 0.5)",
		(x -> (x - 0.5)^2) => "(ξ - 0.5)²",
		(x -> (x - 0.5)^3) => "(ξ - 0.5)³",
		(x -> (x - 0.5)^4) => "(ξ - 0.5)⁴",
		exp => "exp(ξ)", 
		(x -> x * sinpi(x)) => "ξ * sin(ξ * π)",
	])
	numint_nqpoints_slider = @bind numint_nqpoints Slider(1:8; default = 1, show_value = true)
	md"""
	### Numerical integration
	In order to establish the linear equation system above, we have to integrate over the domain [-1, 1]. While we can do this by hand in many cases, but in general this will be a lot of work and in many more advanced cases not possible. Therefore, we typically employ numerical integration in FEM, specifically Gauss Quadrature. To do that, we approximate an integral as a sum over $N_\mathrm{qp}$ quadrature points,
	```math
	\begin{align}
	    I_\mathrm{true} = \int_{-1}^{1} h(\xi)\ \mathrm{d}\xi\ \approx\ \sum_{q = 1}^{N_\mathrm{qp}} h(\xi_q) w_q = I_\mathrm{approx}
	\end{align}
	```
	where the points, $\xi_q$, and weights, $w_q$, are tabulated in the formula sheet.

	To illustrate this, you can pick which function to integrate and how many quadrature points to use below:

	Select function, ``h(\xi) = `` $(numint_fun_selector)

	Number of points: $(numint_nqpoints_slider)

	In the following plot, the width of each bar is the weight, ``w_q``, and the height ``h(\xi_q)``. The marker points show the coordinates ``(\xi_q, h(\xi_q))``. 
	"""
end

# ╔═╡ 1051d025-50fe-4311-94a7-606c46a9e84d
begin
	function integrate(f::F, nqp::Int) where {F<:Function}
		qr = QuadratureRule{RefLine}(nqp)
		return sum(qr.weights .* f.(first.(qr.points)))
	end
	function update_numint_plot!(data, f::F, nqp::Int) where {F<:Function}
		qr = QuadratureRule{RefLine}(nqp)
		ξ_qp = first.(qr.points)
		Plt.update!(data.barplot; arg1 = ξ_qp, arg2 = f.(ξ_qp), width = qr.weights)
		Plt.update!(data.lineplot; arg2 = f.(data.ξv))
		Plt.update!(data.scatterplot; arg1 = ξ_qp, arg2 = f.(ξ_qp))
		Iapprox = integrate(f, nqp)
		Itrue = integrate(f, 12)
		Iapprox_str = @sprintf("%0.5f", Iapprox)
		Itrue_str = @sprintf("%0.5f", Itrue)
		err_str = @sprintf("%0.3e", Iapprox - Itrue)
		data.title[] = L"\text{error} = I_\mathrm{approx} - I_\mathrm{true} = %$(Iapprox_str) - (%$(Itrue_str)) = %$err_str"
		#data.title[] = L"\text{relative error} = 1 - \frac{I_\mathrm{approx}}{I_\mathrm{true}} = 1 - \frac{%$(Iapprox_str)}{%$(Itrue_str)} = %$err_str"
	end
	update_numint_plot!(data_numint, numint_fun, numint_nqpoints)
	fig_numint
end

# ╔═╡ 889968e5-2713-4f0d-9b7c-b83a9180f2ae
md"""
For *Gauss Quadrature* with a certain number quadrature points, weights and locations have been optimized to integrate an as high order polynomial as possible exact. Notice that the linear function is exact for a single point, whereas already with 2 quadrature points, both the 2nd and 3rd order polynomials are integrated exactly. 3 points are needed for the 4th order polynomial. For the exponential and sinusoidal based functions, we observe that the error rapidly decrease with the number of quadrature points.

Equipped with the ability to integrate any function numerically, we can now calculate the integrals in L1a, to establish ``K_{ij}`` and ``f_i``. Specifically, if we use the `linear_line_reference_shape_values` or `quadratic_line_reference_shape_values` functions in `BasicFEM`, we obtain for a given coordinate ``\xi``, the vector ``\underline{N} = [\hat{N}_1(\xi), \cdots, \hat{N}_{N_\mathrm{s}}(\xi)]``.

So we would like to calculate,
```math
    K_{ij} = \int_{-1}^{1} \hat{N}_i(\xi) \hat{N}_j(\xi)\ \mathrm{d}\xi 
	\quad \text{and} \quad
    f_i = \int_{-1}^{1} \hat{N}_i(\xi) g(\xi)\ \mathrm{d}\xi
```

We can then calculate our matrix `Ke` and vector `fe` as (we assume they have been initialized as zeros)
```
[weights, points] = line_quadrature(nquadpoints);
for q = 1:length(weights)
    xi = points(q);
    w = weights(q);
    N = linear_line_reference_shape_values(xi);
    g_xi = g(xi);
    for i = 1:size(N, 2)
        fe(i) = fe(i) + N(i) * g_xi * w;
        for j = 1:size(N, 2)
            Ke(i, j) = Ke(i, j) + N(i) * N(j) * w;
        end
    end
end
```
However, the nested for-loops are quite slow in `MATLAB`, and by using the vector-operations directly, we can calculate this more efficiently as
```
[weights, points] = line_quadrature(nquadpoints);
for q = 1:length(weights)
    xi = points(q);
    w = weights(q);
    N = linear_line_reference_shape_values(xi);
    g_xi = g(xi);
	fe = fe + N' * g(xi) * w;
	Ke = Ke + N' * N * w;
end
```
which will be much faster. Note that in *compiled* languages (e.g. `c++`, `fortran`, `julia`, etc.), using loops can be equally efficient as these vectorized multiplications. 

These implementations are part of the first week homework tasks. 
"""

# ╔═╡ e134610b-dd57-4f67-bc00-b5207288d381
md"""
## L2: Introduction to MATLAB
See notes on Canvas
"""

# ╔═╡ 4e1b0052-2f20-453e-a259-b3e9f1b62752
begin
	parameteric_elem1d_x1_slider = @bind parameteric_elem1d_x1 Slider(0:0.2:2; default = -1.2, show_value = true)
	parameteric_elem1d_x2_slider = @bind parameteric_elem1d_x2 Slider(2.2:0.2:4; default = 1.2, show_value = true)
	parameteric_elem1d_ξ_slider = @bind parameteric_elem1d_ξ Slider(-1:0.1:1; default = 0.0, show_value = true)
	function setup_parametric_elem1d_figure()
		fig = Plt.Figure(size = (600, 100))
		ax_ref = Plt.Axis(fig[1,1]; xlabel = L"\xi")
		ξ = Plt.Observable([NaN]);
		Plt.lines!(ax_ref, [-1, 1], zeros(2); color = :black, linewidth = 2)
		Plt.scatter!(ax_ref, ξ, zeros(1); color = :red)
		Plt.hideydecorations!(ax_ref)
		ax_glob = Plt.Axis(fig[1,2]; xlabel = L"x")
		Plt.xlims!(ax_glob, 0.0, 4.0)
		x = Plt.Observable([NaN, NaN])
		Plt.lines!(ax_glob, x, zeros(2); color = :black, linewidth = 2)
		Plt.scatter!(ax_glob, x, zeros(2); color = :black)
		x_ξ = Plt.Observable([NaN])
		Plt.scatter!(ax_glob, x_ξ, zeros(1); color = :red)
		Plt.hideydecorations!(ax_glob)
		return fig, (; ξ, x, x_ξ)
	end
	parametric_elem1d_fig, parametric_elem1d_data = setup_parametric_elem1d_figure();
	md"""
	## L3a: Parametric elements (1D)
	In order to allow for arbitrary element shapes and coordinates, while retaining the simple description on a fixed **reference shape**, the concept of *parametric elements* is usually employed in finite element analyses. For such elements, the spatial coordinate, ``x``, is described as a function of the reference coordinate, ``\xi``. Specifically, we use shape functions to describe the spatial coordinate with the nodal coordinates, ``x_\alpha``, as coefficients:
	"""
end

# ╔═╡ 40615160-0e12-4a1c-877a-507d1c7b834b
md"""
	The **reference line**, is defined as the line between ``\xi = -1`` and ``\xi = +1`` (which, not by chance 😉, coincides with the domain we used in Lecture 1). The mapping is illustrated for a 2-noded line element (`LinearLine`) below. You can choose the node coordinates, ``x_1`` and ``x_2``, as well as the local coordinate, ``\xi``, resulting in the global coordinate, ``x``.

	``x_1 = `` $(parameteric_elem1d_x1_slider)

	``x_2 = `` $(parameteric_elem1d_x2_slider)

	``ξ = `` $(parameteric_elem1d_ξ_slider)
	"""

# ╔═╡ abf84c67-87f2-4b94-b5a3-578b854e0711
let x1 = parameteric_elem1d_x1, x2 = parameteric_elem1d_x2, ξ = parameteric_elem1d_ξ
	tostr(x) = @sprintf("%4.2f", x)
	N1 = (1 - ξ) / 2
	N2 = (1 + ξ) / 2
	x = N1 * x1 + N2 * x2
	md"""
	``x = \hat{N}_1(\xi) x_1 + \hat{N}_2(\xi) x_2 = ``$(tostr(N1)) * $(tostr(x1)) + $(tostr(N2)) * $(tostr(x2)) = $(tostr(x))
	"""
end

# ╔═╡ d55d3180-eb0b-4d8b-afb2-4daa3dee7b65
begin
	function update_parametric_elem1d_figure!(data, x1, x2, ξ)
		data.ξ[] = [ξ]
		data.x[] = [x1, x2]
		ip = Lagrange{RefLine,1}()
		N = zeros(2)
		Ferrite.reference_shape_values!(N, ip, Vec{1}((ξ,)))
		data.x_ξ[] = [N[1] * x1 + N[2] * x2]
		return data
	end
	update_parametric_elem1d_figure!(parametric_elem1d_data, parameteric_elem1d_x1, parameteric_elem1d_x2, parameteric_elem1d_ξ)
	parametric_elem1d_fig
end

# ╔═╡ c3d23889-ff92-476d-92d0-726812079777
md"""
### Numerical integration
When solving a FEM problem, we end up with **integrals in the spatial domain**,
```math
    \int_a^b h(x)\ \mathrm{d}x
```
When integrating over a single element, the coordinates of the nodes are then ``a = x(\xi=-1)`` and ``b = x(\xi=1)``. We would then like to transform into an integral on the reference domain, ``\xi\in[-1, 1]``. To do this, we simply perform integration by substitution,
```math
    \int_a^b h(x)\ \mathrm{d}x = \int_{x(\xi=-1)}^{x(\xi=1)} h(x)\ \mathrm{d}x = \int_{-1}^{1} h(x(\xi)) x'(\xi)\ \mathrm{d}\xi 
```
The derivative, ``x'(\xi)``, is denoted the jacobian, ``J(\xi)``, and can be calculated as 
```math
    J(\xi) := x'(\xi) = \frac{\mathrm{d}x}{\mathrm{d}\xi} = \sum_{\alpha=1}^{N_\mathrm{nodes}} \frac{\mathrm{d}\hat{N}_\alpha}{\mathrm{d}\xi} x_\alpha
```

Equipped with this result, we can therefore **numerically integrate** a function, ``h(x)``, on the physical domain, ``x \in [a, b]``, by using the quadrature rule on the **reference domain**, ``\xi \in [-1, 1]``, as
```math
    \int_a^b h(x)\ \mathrm{d}x = \int_{-1}^{1} h(x(\xi)) J(\xi)\ \mathrm{d}x \approx \sum_{q = 1}^{N_\mathrm{q}} h(x(\xi_q)) J(\xi_q) w_q
```
where we use the previously introduced equation to calculate the spatial coordinate, ``x_q = x(\xi_q)``, and the jacobian, ``J(\xi_q)``, for each quadrature point, ``\xi_q``. Note that in many cases, such when evaluating shape functions and their gradients, the function ``h(x)`` may be directly expressed as a function of the reference coordinate, i.e. ``g(\xi) = h(x(\xi))`` and we can use ``g(\xi_q)`` directly without calculating ``x_q``. 
"""

# ╔═╡ 6fd9e877-66dc-4a59-8fe2-6d6c9c5af473
begin
	function setup_parametric_line_gradient_fig()
		fig = Plt.Figure(size = (700, 500))
		ax_ref = Plt.Axis(fig[1,1]; xlabel = L"\xi", ylabel = L"\hat{N}(\xi)")
		Plt.lines!(ax_ref, [-1, 1], [1, 0]; label = L"\hat{N}_1(\xi) = [1-\xi]/2,\ \partial\hat{N}_1/\partial\xi = -1/2")
		Plt.lines!(ax_ref, [-1, 1], [0, 1]; label = L"\hat{N}_2(\xi) = [1+\xi]/2,\ \partial\hat{N}_2/\partial\xi = +1/2")
		Plt.ylims!(ax_ref, -0.1, 1.8)
		Plt.xlims!(ax_ref, -1.1, 1.1)
		Plt.Legend(fig[1,2], ax_ref; tellheight=false, tellwidth=false)
		#Plt.axislegend(ax_ref; position = :ct)
		ax_glob = Plt.Axis(fig[2,1:2]; xlabel = L"x", ylabel = L"N(x)")
		Plt.xlims!(ax_glob, -0.1, 4.4)
		Plt.ylims!(ax_glob, -0.1, 1.8)
		N1_line = Plt.lines!(ax_glob, [0.0, 1.0], [1.0, 0.0]; label = L"N_1(x),\ \partial N_1/\partial x = -1/L_e")
		N2_line = Plt.lines!(ax_glob, [0.0, 1.0], [0.0, 1.0]; label = L"N_2(x),\ \partial N_2/\partial x = 1/L_e")
		Plt.axislegend(ax_glob; position = :ct)
		return fig, (;N1_line, N2_line)
	end
	parametric_line_gradient_fig, parametric_line_gradient_data = setup_parametric_line_gradient_fig();
	parametric_line_length_slider = @bind parametric_line_length Slider(0.1:0.1:4.0; default = 2, show_value = true)
md"""
### Mapping of gradients
When solving actual Partial Differential Equations (PDEs), we will also need to get the gradient of a function, ``f(x)``, that we approximate using shape functions, ``g(x) \approx f(x; \underline{a}) = \sum_{i=1}^{N_\mathrm{s}} \hat{N}_i(\xi(x)) a_i``. Specifically, we want to calculate
```math
    \frac{\partial g}{\partial x} \approx \frac{\partial f}{\partial x} = \sum_{i = 1}^{N_\mathrm{s}} \frac{\partial \hat{N}_i(\xi(x))}{\partial x} a_i
```
Since our shape functions, ``\hat{N}(\xi)``, are described as a function of the reference coordinate, ``\xi``, we now need the inverse mapping: ``\xi = \xi(x)`` if we want to calculate the physical coordinate dependent shape function, ``N_i(x) = \hat{N}_i(\xi(x))`` (Notice that we remove the hat, ``\hat{\bullet}``, when denoting the shape function described by the physical coordinates). However, we don't have an explicit function for this inverse mapping, and instead we use
```math
    1 = \frac{\partial x}{\partial x} = \frac{\partial x}{\partial \xi}\frac{\partial \xi}{\partial x} \Rightarrow \frac{\partial \xi}{\partial x} = \left[ \frac{\partial x}{\partial \xi} \right]^{-1} = J^{-1} 
```
to obtain the derivative, ``\partial \xi/\partial x``. The gradient of a shape function can be calculated as,
```math
    \frac{\partial N_i}{\partial x} = \frac{\partial \hat{N}_i}{\partial \xi} \frac{\partial \xi}{\partial x} = \frac{\partial \hat{N}_i}{\partial \xi} J^{-1}
```
This means that we can now calculate the reference shape gradients, ``\partial N_i/ \partial\xi``, and use those to obtain the spatial gradients when required.

For a **linear line element**, the jacobian, ``J``, is given by
```math
J = \frac{\partial \hat{N}_1}{\partial\xi}x_1 + \frac{\partial \hat{N}_2}{\partial\xi}x_2 = \frac{-1}{2}x_1 + \frac{1}{2}x_2 = \frac{1}{2}[x_2 - x_1] = \frac{L_e}{2}
```
where ``L_e = x_2 - x_1`` is the element length. Hence, we get with our mapping that 
```math
\frac{\partial N_i}{\partial x} = \frac{2}{L_e}\frac{\partial \hat{N}_i}{\partial \xi}
```

By changing the element length, we then naturally see that for a short element, the gradient of the shape functions in the physical space will be high, whereas for a longer element, the gradient will be lower. And this is described by the changing jacobian.

Element length, ``L_e =``$(parametric_line_length_slider)
"""
end

# ╔═╡ da8d686b-ccdf-49b0-9fa0-20b2c5b339c8
md"""
``J = ``$(parametric_line_length/2)

``J^{-1} =`` $(round(2/parametric_line_length; digits = 3))

``\partial N_2/\partial x = `` $(round(1/parametric_line_length; digits = 3))
"""

# ╔═╡ 8b268765-0e20-49fc-b5f8-e12f6f6136b5
begin
	Plt.update!(parametric_line_gradient_data.N1_line; arg1 = [0.0, parametric_line_length])
	Plt.update!(parametric_line_gradient_data.N2_line; arg1 = [0.0, parametric_line_length])
	parametric_line_gradient_fig
end

# ╔═╡ 9e31346f-75ce-4d57-9a3a-644e14d5b7bf
md"""
## L3b: The weak form
To introduce the weak form, we will revisit the problem of approximating a function, ``g(x)``. Since we now know how to work in the physical domain, ``x``, we will do that instead of the reference domain, ``\xi``, but this doesn't change the steps.

In this context, we start without introducing the finite element approximation, but rather stating the so-called **strong form** of what we are looking for, namely that we are looking for ``f(x)`` such that
```math
f(x) = g(x),\quad x \in [a, b]
```

Now, we will perform **the key trick** that we will use many times in this course: We take our strong form, and perform the following steps
1) Multiply by an arbitrary test function, ``\delta f``
2) Integrate over the domain

This leads us to the equation,
```math
\int_{a}^{b} \delta f(x) \left[f(x) - g(x)\right]\ \mathrm{d}x = 0, \quad \forall \delta f
```
which is equivalent to the strong form, as ``\delta f(x)`` can be any function. 

!!! warning "The following demonstrates a wrong equation"
	Notice that **without** multiplying with the arbitrary function, ``\delta f``, integration does not lead to an equivalent formulation, i.e.,
	```math
	\int_{a}^{b} \left[f(x) - g(x)\right]\ \mathrm{d}x = 0
	```
	only states that the the average of ``f(x)`` is equal to the average of ``g(x)``. For example, if ``g(x) = x``, then ``f(x) = [a+b]/2`` is one (out of many) possible solutions

"""

# ╔═╡ a1df8715-af64-4701-9b1e-b5ae42e876e4
md"""
However, the computer is not good at finding with mathematical functions, it is much better at working with numerical values. Therefore, we will have to approximate both ``f(x)`` and ``\delta f(x)``. In this course, we will of course use the finite element shape functions to do this approximation, specifically we do
```math
\begin{align}
f(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) a_i = N_i(x) a_i \\
\delta f(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) c_i = N_i(x) c_i
\end{align}
```
Our problem now becomes to find the coefficients ``a_i``, which ensures that the equation,
```math
\int_{a}^{b} c_i N_i(x) \left[N_j(x) a_j - g(x)\right]\ \mathrm{d}x = 0
```
holds for any values of the coefficients, ``c_i``. Since ``c_i`` doesn't depend on the coordinates, we can move these coefficients outside the integral to obtain
```math
\begin{align}
c_i &r_i = 0 \\
&r_i = \int_{a}^{b} N_i(x) \left[N_j(x) a_j - g(x)\right]\ \mathrm{d}x = 0
\end{align}
```
And since we want this to hold for any values of ``c_i``, it has to hold for ``c_1 = 1`` and all others being zero. Similarly for ``c_2 = 1``, and all others zero, and so on. Hence, the only way this can hold true for any value of ``c_i``, is that ``r_i = 0``, and we obtain
```math
\underbrace{\int_{a}^{b} N_i(x) N_j(x)\ \mathrm{d}x}_{K_{ij}}\ a_j = \underbrace{\int_{a}^{b} N_i(x) g(x)\ \mathrm{d}x}_{f_i}
```
Which is exactly the same we obtained by inserting the approximation directly from the beginning as we did in Lecture 1.
"""

# ╔═╡ 082d93db-9088-4a41-bfc4-d1578622a76c
begin 
md"""
## L3b: The weak form

### Alternative approach using directional derivative:

To introduce a very abstract concept in the finite element analysis, we will start by revisiting the problem of approximating (or fitting) a function. But instead of using the finite element approximation, we will instead consider that we want to find a function, ``f(x)``, which minimizes the L2 error wrt. the known function, ``g(x)``, i.e.,
```math
e = \int_{a}^{b} \frac{1}{2}\left[f(x) - g(x)\right]^2\ \mathrm{d}x
```
However, now the unknown is not numerical coefficients, but a mathematical function. So the question now comes - how can we find the stationary point for this problem? 

The answer is that we look for the function, to which if we add any possible variation, will not decrease the error. So let us consider the new function ``f^\ast(x)``, written as
```math
f^\ast(x)(x) = f(x) + \epsilon \delta f(x)
```
If we now take the derivative of e wrt the amount ``\epsilon``, we add of ``\delta f(x)``, to ``f``, in the limit that ``\epsilon \rightarrow 0``, we get the so-called *directional derivative*,
```math
\lim_{\epsilon \rightarrow 0} \frac{\partial e}{\partial \epsilon} = \int_{a}^{b} \delta f(x) \left[f(x) - g(x)\right]\ \mathrm{d}x
```
As ``e`` should not decrease for any function ``\delta f``, this derivative must be zero. If it is positive, then we could add ``-\delta f`` and the error would decrease, and if it is negative we just add ``\delta f`` and it would decrease as well. Therefore, we require this derivative to be zero for any possible function ``\delta f(x)``. This leads us to the equation,
```math
\int_{a}^{b} \delta f(x) \left[f(x) - g(x)\right]\ \mathrm{d}x = 0
```

**Yes, this is abstract:** We now have the function we are looking for, ``f(x)``, and the arbitrary *variation* of that function, ``\delta f(x)``. The key point, is that if this equation holds for any possible function ``\delta f(x)``, there are no way to reduce the error, ``e``, by modifying the function ``f(x)``. Hence, a solution, ``f(x)``, to the above equation is the best possible function we can find! This is the **weak formulation** of the problem, find ``f(x) = g(x)``, and will yield the equivalent solution. For the problem of "approximating" a function, this seems strange because the initial function is known. However, note that in general, the function ``g(x)`` could be unknown, for example the result of an experiment or a simulation. Due to the *variation* ``\delta f(x)``, this form is sometimes called the *variational form* as well. 

However, the computer is not good at finding with mathematical functions, it is much better at working with numerical values. Therefore, we will have to approximate both ``f(x)`` and ``\delta f(x)``. In this course, we will of course use the finite element shape functions to do this approximation, specifically we do
```math
\begin{align}
f(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) a_i = N_i(x) a_i \\
\delta f(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) c_i = N_i(x) c_i
\end{align}
```
Our problem now becomes to find the coefficients ``a_i``, which ensures that the equation,
```math
\int_{a}^{b} c_i N_i(x) \left[N_j(x) a_j - g(x)\right]\ \mathrm{d}x = 0
```
holds for any values of the coefficients, ``c_i``. Since ``c_i`` doesn't depend on the coordinates, we can move these coefficients outside the integral to obtain
```math
\begin{align}
c_i &r_i = 0 \\
&r_i = \int_{a}^{b} N_i(x) \left[N_j(x) a_j - g(x)\right]\ \mathrm{d}x = 0
\end{align}
```
And since we want this to hold for any values of ``c_i``, it has to hold for ``c_1 = 1`` and all others being zero. Similarly for ``c_2 = 1``, and all others zero, and so on. Hence, the only way this can hold true for any value of ``c_i``, is that ``r_i = 0``, and we obtain
```math
\underbrace{\int_{a}^{b} N_i(x) N_j(x)\ \mathrm{d}x}_{K_{ij}}\ a_j = \underbrace{\int_{a}^{b} N_i(x) g(x)\ \mathrm{d}x}_{f_i}
```
Which is exactly the same we obtained by inserting the approximation directly from the beginning as we did in Lecture 1.
"""
nothing
end

# ╔═╡ dc96f2fd-8afc-40b6-afb0-0a613023f7b7
md"""
## L4: The heat equation (1D)
In heat equation problems, our primary unknown to solve for is the temperature field, ``T(x)``. Our derivation will start with the *1st law of thermodynamics*,

**Energy cannot be created or destroyed**


followed by using Fourier's law for the heat flux,
```math
q = -k \frac{\partial T}{\partial x}
```
Fourier's law states that the heat flux, ``q``, (energy / (time and surface area)) is proportional (with material constant ``k``) to the negative temperature gradient. This implies that energy goes from hot areas to cold areas, as we can observe e.g. when cooking.

### Strong form
To derive the strong form, we need to put the 1st law of thermodynamics into a formula. To do this, we must first define the following quantities

* Internal energy, ``e``: The stored internal energy per volume depends on the temperature, and in the linear case we have that ``e = \rho c_\mathrm{p} T``, where ``\rho`` is the density (mass/volume) and ``c_\mathrm{p}`` the heat capacity (energy / (mass and temperature)). 

* Heat flux, ``q``: The flow of thermal energy through the material (energy / (time and area))

* Heat source, ``h``: Externally supplied heat (e.g. microwave heating) inside the body (energy / (time and volume)).

To derive the strong form of the heat equation, we then consider the grey segment with length, ``\Delta x``, inside a bar with constant area, ``A`` and total length, ``L``.
"""

# ╔═╡ 6e4845f2-a984-43f0-9ef6-39f868ed9e20
LocalResource(joinpath(@__DIR__, "heatequation_viapdf.svg"))

# ╔═╡ 8fffc673-54a3-4c5a-bcae-e4af09cab9bf
md"""
The 1st law of thermodynamics then state that the change of internal energy in our segment is equal to the amount of energy going into the segment via the heat flux, ``q``, and the external heat, ``h``. Specifically, we have that 
```math
A \Delta x \dot{e}(x+\Delta x/2) = h(x + \Delta x/2) A \Delta x + q(x)A - q(x + \Delta x) A
```
where ``\dot{e} = \partial e/\partial t``. We approximate the energy and heat source of the segment by a single quadrature point (which is sufficient for very small ``\Delta x`` values).

First, we divide the equation by ``A \Delta x``, 
```math
\dot{e}(x+\Delta x/2) = h(x + \Delta x/2) - \frac{q(x + \Delta x) - q(x)}{\Delta x}
```
Then, we consider the case when ``\Delta x`` goes to zero, i.e.
```math
\dot{e}(x) = h(x) - \lim_{\Delta x \rightarrow 0} \left[\frac{q(x + \Delta x) - q(x)}{\Delta x}\right] = h(x) - \frac{\partial q}{\partial x}
```
In Lecture 9 we will consider the *transient* heat flow, where ``\dot{e} \neq 0``, but for now we will only consider the stationary case, ``\dot{e} = 0``, which occurs after sufficiently long time has passed with constant loading. We then have the following Differential Equation,
```math
\frac{\partial q}{\partial x} = h(x)
```

However, to solve this problem, we need to know how the heat flux depends on the temperature field. This is unique to the material, but using *Fourier's law*, is the standard assumption. Specifically, we then have
```math
q = -k \frac{\partial T}{\partial x}
```
where ``k`` is a material parameter called the heat conductivity. 


Since this is a partial differential equation, we will also need boundary conditions. For now, we will consider the two most common types

1) Dirichlet Boundary Conditions: We know the value of the primary field, in this case the temperature, ``T``, at the boundary
2) Neumann Boundary Conditions: We know the value of the derivative of the primary field, in this case the heat flux, ``q``, out from the boundary.

"""

# ╔═╡ 6e02d8e9-faef-4f90-af0a-fc3f091b55fc
begin
md"""
### Additional discussion on colloqation point methods
The solution we are looking for, ``T(x)``, is a function. So in order to transform this into a problem that the computer can handle, we would like to introduce the shape functions, ``T(x) \approx N_i(x) a_i``. However, if we insert this directly into the strong form, we get 
```math
\frac{\partial}{\partial x}\left[-k \frac{\partial T}{\partial x}\right] = h(x) \Rightarrow
\frac{\partial}{\partial x}\left[-k \frac{\partial N_i}{\partial x}\right]a_i = h(x)
```
Although this looks like we only have one equation, this is actually an equation for every point inside our domain. If we pick a set of points where we solve this equation, we can actually solve this problem. However, there are two problems with this

1) We need the 2nd derivative, ``\partial^2 N_i/\partial x^2``, but when putting together multiple elements, we will not have a continuous derivative across the element boundaries, and then the 2nd derivative goes to infinity
2) We only satisfy the equation point-wise, and cannot guarantee that the solution deviates strongly away from these points. 

However, such methods can still work well, but not with our shape functions. One example is using *Neural Networks* to approximate the temperature, and train them such that the strong form is fullfilled in each evaluated point. These are called *Physics Informed Neural Network*, but will not be discussed in this course. 

Instead, we will obtain a weak form of our problem, aiming to obtain a form similar to the weak form obtained in Lecture 3. To this end, we notice that if we multiply the strong form by an arbitrary function, ``\delta T(x)``, and integrate over the domain, [a, b], this is equivalent to the original strong form:
```math
\int_{a}^{b} \delta T \frac{\partial q}{\partial x}\ \mathrm{d}x = \int_a^b \delta T\ h(x)\ \mathrm{d}x
```
When doing this step, notice two important points:

1) Without multiplying with the *arbitrary* function, ``\delta T``, integration is not equivalent, i.e.
```math
\int_{a}^{b} \frac{\partial q}{\partial x}\ \mathrm{d}x = \int_a^b h(x)\ \mathrm{d}x
```
only states that the strong form is fulfilled on average over the domain.
"""
nothing
end

# ╔═╡ e60e838b-2f53-45d9-b644-615f7171aa0b
md"""
### Weak form
To introduce the weak form, we now do the two steps as above,
1) Multiply with an arbitrary test function, ``\delta T(x)``
2) Integrate over the domain
This yields,
```math
\int_{a}^{b} \delta T \frac{\partial q}{\partial x}\ \mathrm{d}x = \int_a^b \delta T\ h(x)\ \mathrm{d}x
```
Furthermore, we can now perform integration by parts, i.e. we use that 
```math
\frac{\partial}{\partial x}\left[\delta T q\right] = \delta T \frac{\partial q}{\partial x} + \frac{\partial \delta T}{\partial x}q
```
resulting in 
```math
\int_{a}^{b} \delta T \frac{\partial q}{\partial x}\ \mathrm{d}x 
 = \int_{a}^{b} \frac{\partial}{\partial x}\left[\delta T\ q\right]\ \mathrm{d}x - \int_{a}^{b} \frac{\partial \delta T}{\partial x}q\ \mathrm{d}x
= \Big[ \delta T\ q \Big]_a^b - \int_{a}^{b} \frac{\partial \delta T}{\partial x}q\ \mathrm{d}x
```

Putting it all together then yields the weak form of the heat equation in 1D:
```math
- \int_{a}^{b} \frac{\partial \delta T}{\partial x}q\ \mathrm{d}x = \int_a^b \delta T\ h(x)\ \mathrm{d}x - \Big[ \delta T\ q \Big]_a^b, \quad \forall\ \delta T
```

Now we note that if we have Boundary Conditions (BCs) of Neumann type (``q`` is known), we can insert these directly. When we have Dirichlet BCs, the corresponding fluxes will be unknown, and we will deal with this when solving the FE problem. A complete formulation of the FE problem must then state both the weak form including the given Neumann BCs, as well as the Dirichlet BCs as separate equations. 

!!! note "Example"
	If the temperature, ``T_b``, is known at ``x = b``, and the flux, ``q_a``, is known at ``x = a``, we have the complete weak form:

	Find ``T(x)`` such that
	```math
	\begin{align}
	- \int_{a}^{b} \frac{\partial \delta T}{\partial x}q\ \mathrm{d}x &= \int_a^b \delta T\ h(x)\ \mathrm{d}x - \delta T(b)\ q(b) + \delta T(a) q_a, \quad \forall\ \delta T \\
	T(b) &= T_b
	\end{align}
	```	

"""

# ╔═╡ 7ee86550-b9cd-4f6e-b039-c59d9c50f563
md"""
### FE form
To obtain the Finite Element form, we simply use the same technique as earlier, we want to approximate both the test, ``\delta T(x)``, and trial, ``T(x)``, functions using our shape functions:
```math
\begin{align}
\delta T(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) c_i = N_i(x) c_i \\
T(x) &\approx \sum_{i=1}^{N_\mathrm{s}} N_i(x) a_i = N_i(x) a_i
\end{align}
```
Inserting this into the weak form, we obtain
```math
- \int_{a}^{b} \frac{\partial N_i c_i}{\partial x}q\ \mathrm{d}x = \int_a^b N_i c_i\ h(x)\ \mathrm{d}x - \Big[ N_i c_i\ q \Big]_a^b, \quad \forall\ c_i
```
and using that ``c_i`` are arbitrary and not depending on the coordinate, we have
```math
- \int_{a}^{b} \frac{\partial N_i}{\partial x}q\ \mathrm{d}x = \int_a^b N_i\ h(x)\ \mathrm{d}x - \Big[ N_i\ q \Big]_a^b
```
using the same argument as before that ``c_i r_i = 0 \Rightarrow r_i = 0`` if ``c_i`` are arbitrary.

Now you might wonder - we never used ``T(x) \approx N_j(x) a_j``? The reason we cannot see it, is because it enters into the flux ``q`` via Fourier's law, ``q = -k \partial T/\partial x``. Inserting this yields,
```math
\underbrace{\int_{a}^{b} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x}_{K_{ij}}\ a_j = \underbrace{\int_a^b N_i\ h(x)\ \mathrm{d}x + N_i(a) q(a) - N_i(b) q(b)}_{f_i}
```
However, note that we **do not insert the constitutive law in the boundary term!** The reason is that we later want to be able to postprocess to calculate the boundary flux, and inserting the constitutive law at this stage would make that more difficult.
"""

# ╔═╡ 177135bc-71b8-4298-9c21-8c8a11c0c753
md"""
### Solving the FE problem
Before being able to solve the problem, we need to see how the boundary conditions affect it. Specifically, the unknown values, ``a_i``, cannot take any value, they must respect the Dirichlet boundary conditions. If we evaluate the temperature approximation at ``x = a``, we have ``T(a) = N_1(a) a_1 + N_2(a) a_2 + \cdots N_{N_s}(a) a_{N_s}``. However, as mentioned earlier, at the node corresponding to the shape function number, in this case node nr 1 has coordinate ``x = a``, we have that ``N_1(a) = 1``, and ``N_i(a) = 0`` for ``i>1``. Hence, ``T(a) = a_1``, and we must therefore constrain ``a_1 = T_a`` before solving the equation system. We will do this by defining the constrained dofs, `cdofs` (in this case `cdofs = [1]`), and their corresponding known values, `ac` (in this case `ac = [Tₐ]`). Correspondingly, the remaining dofs are free: `fdofs` (in this case, `fdofs = [2, 3, ..., Nₛ]`) with unknown values `af`.  

For the sake of an example with a quadratic line element, we have ``N_\mathrm{s} = 3``, and then `fdofs = [2,3]`. We then want to *partition* our equation system, such that we have
```math
\underline{\underline{K}}\ \underline{a} = \underline{f} \Rightarrow \begin{bmatrix}
\underline{\underline{K}}_{ff} & \underline{\underline{K}}_{fc} \\
\underline{\underline{K}}_{cf} & \underline{\underline{K}}_{cc}
\end{bmatrix}
\begin{bmatrix} \underline{a}_f \\ \underline{a}_c \end{bmatrix} = \begin{bmatrix} \underline{f}_f \\ \underline{f}_c \end{bmatrix}
```
where we in `MATLAB` can obtain these submatrices and vectors as 
```
K_ff = K(fdofs, fdofs); K_fc = K(fdofs, cdofs);
K_cf = K(cdofs, fdofs); K_cc = K(cdofs, cdofs);
a_f = a(fdofs); a_c = a(cdofs);
f_f = f(fdofs); f_c = f(cdofs);
```
"""

# ╔═╡ a004f4bc-1a6a-44f5-810b-58d446f813f4
md"""
### Reaction flux on Dirichlet boundaries
At Dirichlet boundaries, the flux ``q`` is unknown, causing the integral 
```math
f_i = \int_a^b N_i\ h(x)\ \mathrm{d}x + N_i(a) q(a) - N_i(b) q(b)
```
to be unknown. However, for a node-coordinate ``x_j``, we have
```math
N_i(x_j) = \left\lbrace \begin{matrix} 1 & i = j \\ 0 & i \neq j \end{matrix}\right.
```
For the boundary with coordinate ``b``, assume we have node number ``n_b``, such that we have
```math
f_{n_b} = \int_a^b N_{n_b}\ h(x)\ \mathrm{d}x + \underbrace{N_{n_b}(a) q(a)}_{N_{n_b}(a) = 0} - \underbrace{N_{n_b}(b) q(b)}_{N_{n_b}(b) = 1} 
= \int_a^b N_{n_b}\ h(x)\ \mathrm{d}x - q(b)
```
The flux across the boundary ``b``, is then given as 
```math
q(b) = \int_a^b N_{n_b}\ h(x)\ \mathrm{d}x - f_{n_b}
```

In practice, in `MATLAB`, we solve this as
```
[K, f] = calculate_matrix_and_vector(...) % Includes heat source contribution to f
a = zeros(ndofs, 1);
a(cdofs) = ac; % Set constrained values
a(fdofs) = K(fdofs, fdofs) \ (f(fdofs) - K(fdofs, cdofs) * ac); % Solve eq. system
r = zeros(ndofs, 1);
r(cdofs) = f(cdofs) - K(cdofs, :) * a; % Calculate reaction fluxes
```
Where the reason for putting the reaction forces in `r(cdofs)` with preallocated `r = zeros(ndofs, 1)`, is to be able to get index reaction by the global dof number. 
"""

# ╔═╡ 3a99bad7-8447-435d-ae76-392361ca656d
md"""
## L5: Multiple elements - mesh
The true strength of the finite element method comes from combining multiple elements into a mesh (or sometimes called grid). It allows us to approximate the sought function by functions defined on patches of our domain, we will call these patches for *elements*. Consider the heat equation problem above, then we need to calculate 
```math
\begin{align}
K_{ij} &= \int_{a}^{b} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x \\
f_i &= \int_a^b N_i\ h(x)\ \mathrm{d}x + N_i(a) q(a) - N_i(b) q(b)
\end{align}
```
Focusing on the matrix ``K_{ij}``, we note that if we split our domain into ``N_e`` elements, with start and endpoints ``a_e`` and ``b_e`` (where ``e`` is the element number), we can equally well evaluate this integral as a sum of the contributions for each element:
```math
K_{ij} = \int_{a}^{b} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x = \sum_{e = 1}^{N_e} \int_{a_e}^{b_e} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x
```
Before discussing this further, let's discuss the mesh data structure, considering two linear line elements as an example:
"""

# ╔═╡ 9b75b720-0924-4933-9e60-3a0662dfeea7
LocalResource(joinpath(@__DIR__, "mesh_1d_viapdf.svg"))

# ╔═╡ de47e61e-f6db-42b2-96a7-8cd8f31d8825
md"""
To define this mesh, we start by defining the **nodes**. All we need to know about each node is its coordinates, so we define
```
node_coordinates = [x₁, x₂, x₃]
```
Next, we define the **elements**. Here we need to know which nodes belong to each element:
```
e1_nodes = [1; 2]
e2_nodes = [2; 3]
element_nodes = [e1_nodes, e2_nodes]
```

When working with finite element codes, we then have to differentiate between the *local* (element) level and the *global* (system) level. On the local level, we have the node coordinates $x_1^{(e)}$ & $x_2^{(e)}$, and the corresponding shape functions, $N_1^{(e)}$ & $N_2^{(e)}$, where $e$ denotes the element number. 
To get the node coordinates of the element nr `e`, we then do
```
enods = element_nodes(:, e)
ecoords = node_coordinates(:, enods)
```
This is required to be able to evaluate the contributions from the integrals we previously had from the element `e`.

On the global level, we have the node coordinates $x_1$, $x_2$, and $x_3$, with the associated shape functions $N_1(x)$, $N_2(x)$, and $N_3(x)$. Here, we notice that the global shape function $N_2(x)$ is defined by
```math
N_2(x) = \left\lbrace 
\begin{matrix} 
N_2^{(1)}(x) & \text{if }x \in e_1 \\ 
N_1^{(2)}(x) & \text{if }x \in e_2 \\
0 & \text{else}
\end{matrix} \right.
```
Consider when we use the function you have written for the 1D heat equation: `linear_line_heat_element`. This function calculates the *local* stiffness `Ke` and load vector `fe`: 
```math
\begin{align}
K^{(e)}_{ij} &= \int_{x_1^{(e)}}^{x_2^{(e)}} \frac{\partial N_i^{(e)}}{\partial x}\ k\ \frac{\partial N_j^{(e)}}{\partial x}\ \mathrm{d}x\\
f^{(e)}_{i} &= \int_{x_1^{(e)}}^{x_2^{(e)}} N_i^{(e)}\ h(x)\ \mathrm{d}x
\end{align}
```
These contributions are indexed by the *local* indices (in this case $i$ and $j$ are [1, 2]). Therefore, we must know where to put the contributions in the global matrix and vector, ``K_{ij}`` and ``f_i`` (which have indices $i$ and $j$ in [1, 2, 3]). For a calculation of the element `e`, we then do the following
```
[Ke, fe] = linear_line_heat_element(k, b, ecoords, nquadpoints);
edofs = enods;
K(edofs, edofs) = K(edofs, edofs) + Ke;
f(edofs) = f(edofs) + fe;
```

For scalar problems, we will choose to define the shape function number, or Degree of Freedom (DoF) number, to be equal to the node number, i.e. `edofs = enods`. For vector-valued problems, this is not possible, but we will get back to how to solve this later in the course.

The process of calculating and adding the local contributions to the global matrix and vector is called **assembly**, and the complete code for this given a mesh (`element_nodes` and `node_coordinates`), the material parameter `k`, the volumetric heat source `b`, and the number of quadrature points, `nquadpoints`, is
```
nel = size(element_nodes, 2);        % Number of elements
ndofs = size(node_coordinates, 2);   % Number of nodes
K = spalloc(ndofs, ndofs, 3 * ndofs);% Allocate global matrix
f = zeros(ndofs, 1);				 % Allocate global vector
for e = 1:nel
	enods = element_nodes(:, e);     % Extract node numbers for given element
	ecoords = node_coordinates(:, enods); % Extract node coordinates for elem
    % Calculate the local matrix and vector for the given element:
	[Ke, fe] = linear_line_heat_element(k, b, ecoords, nquadpoints);
	% Add the local contributions to the right locations in the global 
    % matrix and vector:
    edofs = enods;
	K(edofs, edofs) = K(edofs, edofs) + Ke;
	f(edofs) = f(edofs) + fe;
end
```
"""

# ╔═╡ aec4b7af-d986-4d33-9f96-dd6c749df0d3
md"""
### Boundary conditions
When applying boundary conditions, we have to work with the *global* indices. In the mesh structure, we also get `node_sets` and `facet_sets`. We will not use the latter in 1d, but will get back to that in 2D. Here we will only use the `node_sets`. We can generate a mesh with the structure above by calling
```
[element_nodes, node_coordinates, facet_sets, node_sets] = ...
	generate_mesh("LinearLine", 2);
```
In this case, we have two node sets:
* `node_sets{"left"}` (in this case: [1])
* `node_sets{"right"}` (in this case: [3])
which contain the left and right node number.

#### Neumann BC
Even in the global numbering, we always have that the shape function corresponding to a specific node is 1 at that node, and all other shape functions are zero. Therefore, if we have Neumann boundary conditions, (known flux), at the right boundary, we want to add the contribution to the global load vector, specifically the contribution
```math
- N_i(b) q(b)
```
where `b` is the coordinate of the end point, i.e. the right node. Hence, we have that
```math
N_i(b) = \left\lbrace \begin{matrix} 0 & i \neq 3 \\ 1 & i = 3 \end{matrix} \right.
```
such that we can simply add this contribution as,
```
right_node = node_sets{"right"};
f(right_node) = f(right_node) - qb
```
where `qb` is our known boundary flux. 

#### Dirichlet BC
Utilizing the same property, that ``N_i(x_j)`` is 1 if ``i = j``, and 0 if ``i \neq j``, we have that the temperature on the left side, ``T(a) = a_1``, such that we can enforce the dirichlet boundary conditions and solve the equation system via partitioning as 
```
cdofs = node_sets{"left"}
ac = [12] % Set the temperature on the left side to 12 degrees
fdofs = setdiff((1:ndofs)', cdofs) % All other dofs
a = zeros(ndofs, 1)				% Solution vector
a(cdofs) = ac; 					% Set prescribed values
% Solve the equation system
a(fdofs) = K(fdofs, fdofs)\(f(fdofs) - K(fdofs, cdofs) * ac);
```
Now we have solved our first finite element problem with multiple elements 🎉

And the best part: The code works equally well for any number of elements: just change the number of elements passed to the `generate_mesh` function call! 
"""

# ╔═╡ 0363492e-8ecf-4f6c-ab50-2b0612572d39
md"""
## L6: The heat equation (2D)
"""

# ╔═╡ 080f125a-8003-41e7-b50d-05a4de587301
md"""
### Divergence
Consider a vector field, 
```math
\underline{q}(\underline{x}) = \begin{bmatrix} q_1(\underline{x}) \\ q_2(\underline{x}) \end{bmatrix}, \quad \underline{x} = \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
```
in a 2D body with thickness ``t``, containing a the blue rectangle with volume ``V = \Delta x_1\ \Delta x_2\ t``:
"""

# ╔═╡ 58c1a7d8-03fb-409d-8589-46d2fa597f82
LocalResource(joinpath(@__DIR__, "divergence_body_viapdf.svg"))

# ╔═╡ 8786dc4e-f7f4-4ccc-85f5-5302988c696c
md"""
The definition of the divergence of ``\underline{q}`` is
```math
\mathrm{div}(\underline{q}) := \lim_{V\rightarrow 0}\left[ \frac{1}{V} \int_{\Gamma} \underline{q} \cdot \underline{n}\ \mathrm{d}\Gamma \right]
```
If we consider the surface integral, we have
```math
\begin{align}
\frac{1}{t}\int_{\Gamma} \underline{q} \cdot \underline{n}\ \mathrm{d}\Gamma 
=& \int_{x_1}^{x_1 + \Delta x_1} \underline{q}(\tilde{x}, x_2) \cdot[-\underline{e}_2]\ \mathrm{d}\tilde{x}
+ \int_{x_1}^{x_1 + \Delta x_1} \underline{q}(\tilde{x}, x_2 + \Delta x_2) \cdot \underline{e}_2\ \mathrm{d}\tilde{x} \\
+& \int_{x_2}^{x_2 + \Delta x_2} \underline{q}(x_1, \tilde{x}) \cdot [-\underline{e}_1]\ \mathrm{d}\tilde{x}
+ \int_{x_2}^{x_2 + \Delta x_2} \underline{q}(x_1 + \Delta x_1, \tilde{x}) \cdot \underline{e}_1\ \mathrm{d}\tilde{x} \\
=&-\int_{x_1}^{x_1 + \Delta x_1} q_2(\tilde{x}, x_2) \ \mathrm{d}\tilde{x}
+ \int_{x_1}^{x_1 + \Delta x_1} q_2(\tilde{x}, x_2 + \Delta x_2)\ \mathrm{d}\tilde{x} \\
&- \int_{x_2}^{x_2 + \Delta x_2} q_1(x_1, \tilde{x})\ \mathrm{d}\tilde{x}
+ \int_{x_2}^{x_2 + \Delta x_2} q_1(x_1 + \Delta x_1, \tilde{x})\ \mathrm{d}\tilde{x} \\
\end{align}
```
For a continuous field, ``\underline{q}(\underline{x})``, and ``\Delta x_1`` and ``\Delta x_2`` small, we can use a single quadrature point to calculate the integrals (this can be shown by using e.g. a taylor series expansion), such that we get,
```math
\begin{align}
\frac{1}{t}\int_{\Gamma} \underline{q} \cdot \underline{n}\ \mathrm{d}\Gamma 
&= \Delta x_1 \left[ q_2(x_1 + \Delta x_1/2, x_2 + \Delta x_2) - q_2(x_1 + \Delta x_1/2, x_2)\right] \\
&+ \Delta x_2 \left[ q_1(x_1 + \Delta x_1, x_2 + \Delta x_2/2) - q_1(x_1, x_2 + \Delta x_2/2)\right]
\end{align}
```
Divide by ``\Delta x_1\ \Delta x_2``, and let ``\Delta x_1 \rightarrow 0`` and ``\Delta x_2 \rightarrow 0``, such that also ``V = \Delta x_1\ \Delta x_2\ t \rightarrow 0``,
```math
\begin{align}
\lim_{V\rightarrow 0} \frac{1}{V}\int_{\Gamma} \underline{q} \cdot \underline{n}\ \mathrm{d}\Gamma 
&= \lim_{\Delta x_2 \rightarrow 0} \frac{q_2(x_1, x_2 + \Delta x_2) - q_2(x_1, x_2)}{\Delta x_2}\\ 
&+ 
\lim_{\Delta x_1 \rightarrow 0} \frac{q_1(x_1 + \Delta x_1, x_2) - q_1(x_1, x_2)}{\Delta x_1}
\end{align}
```
Upon applying the definition of (partial) derivatives, we then obtain,
"""

# ╔═╡ 5608fe81-9b47-4214-a8ad-0960ca147b27
md"""
### Divergence theorem
We denote the full body from above ``\Omega``, and its boundary ``\Gamma``. Then, we split this body in two parts, ``\Omega_1`` and ``\Omega_2``, with boundaries ``\Gamma_1`` and ``\Gamma_2``, respectively. These boundaries consist of ``\tilde{\Gamma}_i \subset \Gamma`` (i.e. ``\tilde{\Gamma}_i`` is part of ``\Gamma``), and ``\bar{\Gamma}_i`` which is the cut boundary. We would then like to show that
```math
\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_{\Gamma_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\Gamma_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma
```
"""

# ╔═╡ 340f8448-6db6-404f-a7ea-192aa339c286
LocalResource(joinpath(@__DIR__, "divergencethm_split_viapdf.svg"))

# ╔═╡ d9cdd761-bcb3-461a-b473-ab1e6493a80b
md"""
Expanding the integrals, we have 
```math
\begin{align}
\int_{\Gamma_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma &= \int_{\tilde{\Gamma}_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\bar{\Gamma}_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma \\
\int_{\Gamma_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma &= \int_{\tilde{\Gamma}_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\bar{\Gamma}_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma
\end{align}
```
When we note that the cut boundaries, ``\bar{\Gamma}_1`` and ``\bar{\Gamma}_2``, are equal except that the outward pointing normal vector have opposite direction, i.e. ``\bar{\underline{n}}_1(\underline{x}) = -\bar{\underline{n}}_2(\underline{x})``. Consequently, we get
```math
\int_{\bar{\Gamma}_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = -\int_{\bar{\Gamma}_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma
```
such that
```math
\int_{\Gamma_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\Gamma_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_{\tilde{\Gamma}_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\tilde{\Gamma}_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma
```
which is what we wanted to show! 
"""

# ╔═╡ 33618c32-a221-487c-b1ac-0b4beb347c7d
md"""
This argument also holds when dividing into arbitrarily many, ``N_\mathrm{V}``, parts, i.e.
```math
\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma 
= \sum_{i = 1}^{N_V} \left[\int_{\Gamma_i} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma\right]
= \sum_{i = 1}^{N_V} \left[\left[\frac{1}{V_i} \int_{\Gamma_i} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma\right]V_i \right]
```
where ``V_i`` is the volume of part ``i``. If we let the size of each part go to zero, ``V_i \rightarrow 0`` (loosely ``V_i \rightarrow \mathrm{d}\Omega \rightarrow 0``), we then have by using the definition of divergence,
"""

# ╔═╡ ac197737-6f54-4f22-a049-0e268c7ec4db
md"""
Which is called the **divergence theorem**.
"""

# ╔═╡ b1999d4a-3111-4203-ab1c-173905f06d80
md"""
### Strong form
Consider an arbitrary part of a body ``\Omega``. The energy balance is then,
```math
\underbrace{\int_\Omega \dot{e}\ \mathrm{d}\Omega}_{\text{Internal energy change}} = \underbrace{\int_\Omega h\ \mathrm{d}\Omega}_{\text{Volumetric heat supply}} - \underbrace{\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma}_{\text{Boundary heat loss}}
```
Utilizing the divergence theorem on the boundary contribution, we get
```math
\int_\Omega \dot{e}\ \mathrm{d}\Omega = \int_\Omega h\ \mathrm{d}\Omega - \int_\Omega \mathrm{div}(\underline{q})\ \mathrm{d}\Omega
```
This has to hold for any domain, ``\Omega``, so it must hold pointwise, i.e.
```math
\dot{e} = h - \mathrm{div}(\underline{q})
```
This argument is called *localization*, and holds when we only have volume integrals (i.e. no mixed surface and volume integrals). As for the 1D case, we consider the stationary case for now, such that ``\dot{e} = 0``, and we have 
"""

# ╔═╡ 33c4264d-78bb-4c95-9d27-aa340df03808
md"""
As for the 1D case, we will use Fourier's law to model the heat flux, ``q``, given as 
```math
\underline{q} = -\underline{\underline{D}}\ \mathrm{grad}(T)
```
where ``T`` is the temperature, and ``\underline{\underline{D}}`` the conductivity matrix. The gradient, ``\mathrm{grad}(T)``, is
```math
\mathrm{grad}(T) = \begin{bmatrix} \frac{\partial T}{\partial x_1} \\ \frac{\partial T}{\partial x_2} \end{bmatrix}
```
in 2D. The conductivity matrix, ``\underline{\underline{D}}``, is a material parameter, and when the material is *isotropic*, such that the conductivity is the same in all directions, we can write this as ``\underline{\underline{D}} = k\ \underline{\underline{I}}``, where ``\underline{\underline{I}}`` is the identity matrix and ``k`` the scalar heat conductivity. 
"""

# ╔═╡ fccdcef3-3768-47d4-ba0e-ecca7308d356
md"""
### Weak form
To introduce the weak form, we now do the two steps as before,
1) Multiply with an arbitrary test function, ``\delta T(\underline{x})``
2) Integrate over the domain, ``\Omega``
This yields,
```math
\int_\Omega \delta T\ \mathrm{div}(\underline{q})\ \mathrm{d}\Omega = \int_\Omega \delta T\ h\ \mathrm{d}\Omega
```
The integration by parts work similar, but a bit more complicated in 2D, we start by considering
```math
\begin{align}
\mathrm{div}(\delta T\ \underline{q}) &= \frac{\partial}{\partial x_1}\left[\delta T\ q_1\right] + \frac{\partial}{\partial x_2}\left[\delta T\ q_2\right] \\
&=  \frac{\partial \delta T}{\partial x_1}q_1 + \delta T \frac{\partial q_1}{\partial x_1} + \frac{\partial \delta T}{\partial x_2}q_2 + \delta T \frac{\partial q_2}{\partial x_2}\\
&= \mathrm{grad}(\delta T) \cdot \underline{q} + \delta T\ \mathrm{div}(\underline{q})
\end{align}
```
Now we can insert ``\mathrm{div}(\delta T\ \underline{q})`` in the divergence theorem to obtain,
"""

# ╔═╡ 02ebd323-03e2-41cb-903f-cb0bf690cfd5
md"""
Using this in the weak form, we can transform the left hand side,
```math
\int_\Omega \delta T\ \mathrm{div}(\underline{q})\ \mathrm{d}\Omega = \int_\Gamma \delta T\ \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma - \int_\Omega \mathrm{grad}(\delta T)\cdot \underline{q}\ \mathrm{d}\Omega
```
This transforms the integral and introduces a boundary term (similar to integration by parts in 1D), and we get
```math
- \int_\Omega \mathrm{grad}(\delta T)\cdot \underline{q}\ \mathrm{d}\Omega = 
\int_\Omega \delta T\ h\ \mathrm{d}\Omega - 
\int_\Gamma \delta T\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
Here, the boundary flux, ``q_\mathrm{n} := \underline{q}\cdot\underline{n}``, is known for Neumann boundary conditions, and unknown for Dirichlet boundary conditions. 
"""

# ╔═╡ 27087e5b-08fb-4e0e-962c-30999860e238
md"""
### FE form
In the next lectures, we will discuss a lot of the complications when going to 2D (and 3D), but the derivation is very similar to the 1D case: We simply insert the approximations of the test function, ``\delta T(\underline{x})``, and the temperature itself, ``T(\underline{x})``, by using shape functions that are defined for 2D (or 3D) coordinates, ``\underline{x}``,
```math
\begin{align}
\delta T(\underline{x}) &\approx \sum_{i = 1}^{N_\mathrm{s}} N_i(\underline{x}) c_i = N_i(\underline{x}) c_i \\
T(\underline{x}) &\approx \sum_{i = 1}^{N_\mathrm{s}} N_i(\underline{x}) a_i = N_i(\underline{x}) a_i
\end{align}
```
where the coefficients ``c_i`` are arbitrary and ``a_i`` are the unknowns we are looking for (except those known due to Dirichlet BCs). Inserting this into the weak form, we obtain,
```math
- \int_\Omega \mathrm{grad}(N_i(\underline{x}) c_i)\cdot \underline{q}\ \mathrm{d}\Omega = 
\int_\Omega N_i(\underline{x}) c_i\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i(\underline{x}) c_i\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
and noting that ``c_i`` doesn't depend on the coordinate, we can move it outside the (spatial) derivatives and integrals, putting everything on one side of the equal sign, to obtain
```math
c_i \underbrace{\left[-\int_\Omega \mathrm{grad}(N_i(\underline{x}))\cdot \underline{q}\ \mathrm{d}\Omega - 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega + 
\int_\Gamma N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma \right]}_{r_i} = 0
```
Thus we have the form ``c_i r_i = 0``, where ``c_i`` are arbitrary. As before, first we consider ``c_1 = 1``, and all other zero, which gives ``r_1 = 0``. Then ``c_2 = 1`` and all other zero, giving ``r_2 = 0``, and so on, showing that we have ``r_i = 0`` in general, giving the FE form,
```math
-\int_\Omega \mathrm{grad}(N_i(\underline{x}))\cdot \underline{q}\ \mathrm{d}\Omega = 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
Finally, we insert the constitutive relationship - Fourier's law, ``\underline{q} = -\underline{\underline{D}}\ \mathrm{grad}(T(\underline{x})) = - \underline{\underline{D}}\ \mathrm{grad}(N_j(\underline{x}) a_j)``, using that ``a_j`` doesn't depend on the coordinates and can be moved outside (spatial) derivatives and integrals,
```math
\underbrace{\int_\Omega \mathrm{grad}(N_i(\underline{x}))^\mathrm{T} \underline{\underline{D}}\ \mathrm{grad}(N_j(\underline{x})) \ \mathrm{d}\Omega}_{K_{ij}}\ a_j = \underbrace{
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma}_{f_i}
```
"""

# ╔═╡ d9385605-737f-401b-aba9-5d8d7a18a18c
md"""
## L7-8: Finite elements in 2D and 3D

"""

# ╔═╡ a658416a-108e-4f2f-9cb1-324e5f2960f7
md"""
### Shape functions
As for the 1D elements, we define our shape functions on reference elements with a fixed geometry. In 2D, we have 2 different reference elements:
"""

# ╔═╡ 32c1a2ed-ea2a-469b-936f-1a2b4b3eec90
LocalResource(joinpath(@__DIR__, "refshapes_2d.svg"))

# ╔═╡ ef71b38b-ec86-41eb-86f2-ea8087038502
begin
	fig_triangle_shapefuns, data_triangle_shapefuns = let
		function add_axis!(fig, nr)
			Plt.Label(fig[1,nr]; text = L"\hat{N}^{(e)}_{%$nr}", tellwidth = false, fontsize = 28)
			angle = Plt.Observable(1.125π)
			ax = Plt.Axis3(fig[2,nr]; xlabel = L"\xi_1", ylabel = L"\xi_2", zlabel = 	L"\hat{N}(\xi\underbar)", azimuth = angle)
			Plt.limits!(ax, (0, 1), (0, 1), (0, 1))
			ξ1 = [1, 0, 0, 1]
			ξ2 = [0, 1, 0, 0]
			N = zeros(3); N[nr] = 1;
			Plt.lines!(ax, ξ1, ξ2, zeros(4); linewidth = 2, color = :black)
			Plt.scatter!(ax, ξ1, ξ2, zeros(4); color = :black)
			Plt.mesh!(ax, ξ1[1:3], ξ2[1:3], N; alpha = 0.5)
			return angle
		end
		fig = Plt.Figure(size = (800, 300))
		a1 = add_axis!(fig, 1)
		a2 = add_axis!(fig, 2)
		a3 = add_axis!(fig, 3)
		fig, (a1, a2, a3)
	end
	triangle_shapefun_slider = @bind triangle_shapefun_rot Slider(range(0, 2π, 21), default = 1.2π)
	md"""
	Here, it is important to separate between a numbered vector, ``\underline{\xi}_i``, and the component of a vector ``\xi_i``! 
	
	On this reference shapes, we define the vertex coordinates as well as a numbering of each facet. The facets describe the edges - in 3d the facets would be the faces but we don't consider 3d in this course so for us, facet and edge is equivalent. The facets are used when assembling Neumann BCs later.
	
	#### Triangle shape functions
	The linear triangle shape functions are given as 
	```math
	\begin{align}
	    \hat{N}^e_1(\underline{\xi}) = \xi_1, \quad
	    \hat{N}^e_2(\underline{\xi}) = \xi_2, \quad
	    \hat{N}^e_3(\underline{\xi}) = 1 - \xi_1 - \xi_2
	\end{align}
	```
	Rotation: $(triangle_shapefun_slider)
	"""
end

# ╔═╡ 986c99b7-81b8-4ceb-9eee-14e2d0d26b34
let θ = triangle_shapefun_rot
	for i = 1:3
		data_triangle_shapefuns[i][] = θ
	end
	fig_triangle_shapefuns
end

# ╔═╡ d40d0648-ea58-47d1-b2b7-77c30005f50a
md"""
#### Quadrilateral shape functions
Similarly, we define the *bilinear* quadrilateral shape functions as,
```math
\begin{align}
    \hat{N}^e_1(\underline{\xi}) = \frac{[1 - \xi_1][1 - \xi_2]}{4}, \quad
    \hat{N}^e_2(\underline{\xi}) = \frac{[1 + \xi_1][1 - \xi_2]}{4} \\
    \hat{N}^e_3(\underline{\xi}) = \frac{[1 + \xi_1][1 + \xi_2]}{4}, \quad
    \hat{N}^e_4(\underline{\xi}) = \frac{[1 - \xi_1][1 + \xi_2]}{4}
\end{align}
```
#### Multiple elements
When we put these into a mesh, we have the same concept as for the 1D case. For each element, we know its node numbers, e.g. 
```
enods = [3, 4, 7] % Element 5
enods = [4, 8, 7] % Element 6
```
In the mesh below, only elements 5 and 6 share node 4 (Note, it is not possible from the figure to know the element numbers). Hence, the *global* shape number 4 is defined as
```math
N_4(x) = \left\lbrace \begin{matrix} N_2^{(5)}(x) & e = 5 \\ N_1^{(6)}(x) & e = 6 \\ 0 & \text{else}\end{matrix}\right.
```
"""

# ╔═╡ e2f75d06-a29f-4709-a3ca-3dc98cd4f570
begin
	cell_2d3x3_selector = @bind cell_2d3x3_type Select(["Triangle", "Quadrilateral"])
	nodenr_2d3x3_selector = @bind nodenr_2d3x3 Select(collect(1:16); default=4)
	azimuth_3x3_2d_slider = @bind azimuth_3x3_2d Slider(range(0, 2, 21); default=1.2, show_value=true)
md"""
Select cell type: $(cell_2d3x3_selector)

Select shape function from node number: $(nodenr_2d3x3_selector)

Rotation: $(azimuth_3x3_2d_slider) π
"""
end

# ╔═╡ ce1d0468-b97a-4ae1-b93c-28cc10acb7d4
begin
function setup_3x3_2d_grid(::Type{CT}) where {CT<:Union{Triangle, Quadrilateral}}
	grid = generate_grid(CT, (3, 3))
	dh = DofHandler(grid)
	ip = geometric_interpolation(CT)
	add!(dh, :u, ip)
	close!(dh)
	node_to_dof_mapping = zeros(Int, 16)
	node_coordinates = Ferrite.get_node_coordinate.((grid,), 1:16)

	for cell in CellIterator(dh)
		for (node_nr, x_node) in enumerate(node_coordinates)
			@assert length(getcoordinates(cell)) == length(celldofs(cell))
			for (x_cell, dofnr) in zip(getcoordinates(cell), celldofs(cell))
				if x_cell ≈ x_node
					node_to_dof_mapping[node_nr] = dofnr
				end
			end
		end
	end
	
	return dh, node_to_dof_mapping
end
	
makiepoint(x::Vector, y::Vector, z::Vector) = Plt.Point3f.(x, y, z)
makiepoint(x::Number, y::Number, z::Number) = makiepoint([x], [y], [z])
makiepoint(v::Vec{3}) = makiepoint(v...)
makiepoint() = makiepoint(NaN, NaN, NaN)
makiepoint(n::Int) = makiepoint((fill!(zeros(n), NaN) for _ in 1:3)...)

function create_empty_plot_and_observables()
	fig = Plt.Figure()
	azimuth = Plt.Observable(1.275π)
	elevation = Plt.Observable(π/8)
	ax = Plt.Axis3(fig[1,1]; xlabel=L"x_1", ylabel=L"x_2", zlabel=L"N(x)", azimuth, elevation)
	Plt.limits!(ax, (-1, 1), (-1, 1), (0, 1))

	xv = Plt.Observable([-1.0, 1.0])
	yv = Plt.Observable([-1.0, 1.0])
	zv = Plt.Observable([0.0 0.0; 1.0 1.0])
	
	Plt.surface!(ax, xv, yv, zv; colorrange=(-1.0, 1.0))
	
	cell_edges = Plt.Observable(makiepoint())
	nodes = Plt.Observable(makiepoint(16))
	
	Plt.lines!(ax, cell_edges; color=:black, linewidth = 2)
	
	Plt.scatter!(ax, nodes; color=:black, label="Node")
	Plt.text!(ax, nodes; text=string.(collect(1:16)))
	Plt.Legend(fig[1,2], ax)
	observables = Dict(
		"xv"=>xv, "yv"=>yv, "zv"=>zv, 
		"cell_edges"=>cell_edges, "nodes"=>nodes,
		"azimuth" => azimuth, "elevation" => elevation)
	return fig, observables
end

function update_observables!(observables, dh, node_to_dof_mapping, nodenr, azimuth, elevation)
	dofnr = node_to_dof_mapping[nodenr]
	a = zeros(ndofs(dh)); a[dofnr] = 1.0
	numpoints = 100
	xv = collect(range(-1.0, 1.0, numpoints))
	yv = copy(xv)
	points = [Vec(x,y) for y in yv for x in xv]
	ph = PointEvalHandler(dh.grid, points)
	zv = reshape(evaluate_at_points(ph, dh, a), (numpoints, numpoints))
	observables["xv"][] = xv
	observables["yv"][] = yv
	observables["zv"][] = zv
	ncellnodes = Ferrite.nnodes_per_cell(dh.grid, 1) + 2
	cell_edges = makiepoint((ncellnodes) * getncells(dh.grid))
	
	for cell in CellIterator(dh)
		for (i, x_cell) in enumerate(getcoordinates(cell))
			cell_edges[(cellid(cell) - 1) * ncellnodes + i] = makiepoint(x_cell[1], x_cell[2], a[celldofs(cell)[i]])[1]
		end
		cell_edges[cellid(cell) * ncellnodes - 1] = makiepoint(getcoordinates(cell)[1][1], getcoordinates(cell)[1][2], a[celldofs(cell)[1]])[1]
	end
	observables["cell_edges"][] = cell_edges
	
	ncoords = Ferrite.get_node_coordinate.((dh.grid,), 1:16)
	anodes = zeros(length(ncoords)); anodes[nodenr] = 1.0
	observables["nodes"][] = makiepoint(first.(ncoords), last.(ncoords), anodes)
	observables["azimuth"][] = azimuth
	observables["elevation"][] = elevation
	return nothing
end
fig_3x3_2d, obs_3x3_2d = create_empty_plot_and_observables();
end;

# ╔═╡ 76a2f37c-2672-431a-9a33-128681607393
begin
	CT = cell_2d3x3_type == "Triangle" ? Triangle : Quadrilateral
	dh_3x3_2d, node2dof = setup_3x3_2d_grid(CT) 
	update_observables!(obs_3x3_2d, dh_3x3_2d, node2dof, nodenr_2d3x3, azimuth_3x3_2d*π, π/8);
	fig_3x3_2d
end

# ╔═╡ 903164a8-30d5-466f-b399-eb380a895117
md"""
### Parametric elements
As in 1D, we have for the parametric elements that
```math
\hat{N}^{e}_i(\underline{\xi}) = N^{e}_i(\underline{x}(\underline{\xi}))
```
where the parametric mapping is 
```math
\underline{x}(\underline{\xi}) = \sum_{\alpha = 1}^{N_\mathrm{nodes}} \hat{N}^{e}_\alpha(\underline{\xi})\ \underline{x}_\alpha^e
```
where ``\underline{x}_\alpha^e`` is the coordinate of local node number ``\alpha`` in element ``e``.
"""

# ╔═╡ dfe56d4b-13c2-4ad4-ab5e-578f470e4589
md"""
#### Numerical integration
As for 1D, an integration over a 2D reference element, ``\hat{\Omega}``, is approximated by the sum
```math
\int_{\hat{\Omega}} h(\underline{\xi})\ \mathrm{d}\Omega \approx \sum_{q = 1}^{N_\mathrm{qp}} h(\underline{\xi}_q) w_q
```
where ``\underline{\xi}_q`` and ``w_q`` are the ``q``th tabulated integration point and weight. However, for quadrilaterals, we don't tabulate these values, since they can be derived based on the line quadrature. Specifically, we have
```math
\underline{\xi}_q = [\xi_a^{1\mathrm{d}}, \xi_b^{1\mathrm{d}}]^\mathrm{T}, \quad w_q = w_a^{1\mathrm{d}} w_b^{1\mathrm{d}}
```
Given a line quadrature rule with ``N_q^{1\mathrm{d}}`` points, we then use all combinations of ``a \in [1, N_q^{1\mathrm{d}}]`` and ``b \in [1, N_q^{1\mathrm{d}}]`` to create ``N_q = \left[N_q^{1\mathrm{d}}\right]^2`` quadrature points and weights for the quadrilateral. 

The final question for the numerical integration, is then how to modify the weights, ``w_q``, to account for the physical geometry of an element. To this end, we consider how an area element in the reference element, ``\mathrm{d}\hat{\Omega} = \mathrm{d}\xi_1\ \mathrm{d}\xi_2``, transforms as we change to the physical element and the area element ``\mathrm{d}\Omega``,
"""

# ╔═╡ e85161d4-4b28-47d4-a05e-a3c3e34415b5
LocalResource(joinpath(@__DIR__, "numint_map_2d_viapdf.svg"))

# ╔═╡ fdc596a9-bd9e-45bd-a082-0418bb5b5084
md"""
As we see from the physical element, the initial rectangular element ``\mathrm{d}\hat{\Omega}`` is transformed into a parallelogram, with sides ``\underline{v}_1`` and ``\underline{v}_2``. The area is then given by the out of plane component of the cross product, ``\underline{v}_1 \times \underline{v}_2``. Using the Jacobian, ``\underline{\underline{J}}``, we can write the side vectors as
```math
\underline{v}_1 = \begin{bmatrix} J_{11} \\ J_{21} \end{bmatrix} \mathrm{d}\xi_1, \quad 
\underline{v}_2 = \begin{bmatrix} J_{12} \\ J_{22} \end{bmatrix} \mathrm{d}\xi_2
```
And the area becomes
```math
\mathrm{d}\Omega = \underline{e}_3 \cdot \left[\underline{v}_1 \times \underline{v}_2\right] = \left[J_{11}\ J_{22} - J_{21}\ J_{12}\right]\ \mathrm{d}\xi_1\ \mathrm{d}\xi_2 = \mathrm{det}(\underline{\underline{J}})\ \mathrm{d}\hat{\Omega}
```
Based on this, we get the integral transformation
```math
\int_{\Omega} h(\underline{x})\ \mathrm{d}\Omega = \int_{\hat{\Omega}} h(\underline{x}(\underline{\xi})) \mathrm{det}(\underline{\underline{J}})\ \mathrm{d}\hat{\Omega}
```
Similar to before, we can calculate the Jacobian as 
```math
\underline{\underline{J}}(\underline{\xi}) := \frac{\partial \underline{x}}{\partial \underline{\xi}} = \sum_{\alpha = 1}^{N_\mathrm{nodes}^e} \underline{x}_\alpha^e \left[\frac{\partial \hat{N}^e_\alpha}{\partial \underline{\xi}}\right]^\mathrm{T}
```
This becomes a bit more straight-forward using the index notation with summation, where we have 
```math
J_{ij}(\underline{\xi}) := \frac{\partial x_i}{\partial \xi_j} = \sum_{\alpha = 1}^{N_\mathrm{nodes}^e} \left[\underline{x}_\alpha^e\right]_i \frac{\partial \hat{N}^e_\alpha}{\partial \xi_j}
```
In `MATLAB`, given the nodal coordinates as a ``[2, N_\mathrm{nodes}]`` matrix and the derivatives of the shape functions as a ``[2, N_\mathrm{nodes}]`` matrix, we can calculate the sum as a matrix-matrix multiplication. When implementing the `calculate_jacobian` function, validate that if you use
```
coords = [1 0 0; 0 0 0]
dNdxi = [0 0 0; 1 0 0]
```
you get
```
J = [0 1; 0 0]
```
"""

# ╔═╡ 9aec05ec-d7eb-4ef1-a56f-96a4c2371897
md"""
#### Mapping of gradients
As for 1D, we have that
```math
N_i(\underline{x}(\underline{\xi})) = \hat{N}_i(\underline{\xi})
```
However, this implies that if the geometry of the physical element is different from the reference element,
```math
\frac{\partial N_i}{\partial \underline{x}} \neq \frac{\partial \hat{N}_i}{\partial \underline{\xi}}
```
Instead, we have to consider the derivative,
```math
\frac{\partial N_i}{\partial \underline{\xi}} = \frac{\partial N_i}{\partial \underline{x}} \cdot \frac{\partial \underline{x}}{\partial \underline{\xi}}
= \frac{\partial N_i}{\partial \underline{x}} \cdot \underline{\underline{J}},
\quad \text{or in index notation: }\quad 
\frac{\partial N_i}{\partial \xi_j} = \frac{\partial N_i}{\partial x_k} \frac{\partial x_k}{\partial \xi_j} = \frac{\partial N_i}{\partial x_k} J_{kj}
```
Right-multiplying both sides by ``\underline{\underline{J}}^{-1}`` gives
```math
\frac{\partial N_i}{\partial \underline{\xi}} \cdot \underline{\underline{J}}^{-1} = 
\frac{\partial N_i}{\partial \underline{x}},
\quad \text{or in index notation: }\quad
\frac{\partial N_i}{\partial \xi_j} J_{jl}^{-1} = \frac{\partial N_i}{\partial x_k} J_{kj} J_{jl}^{-1} = \frac{\partial N_i}{\partial x_l}
```
or equivalently
```math
\begin{align}
\frac{\partial N_i}{\partial \underline{x}} &= 
\underline{\underline{J}}^{-T} \cdot \frac{\partial N_i}{\partial \underline{\xi}}\\
\frac{\partial N_i}{\partial x_l} &= J_{lj}^{-T} \frac{\partial N_i}{\partial \xi_j} 
\end{align}
```

In `MATLAB`, we will store the gradients as (assuming three shape functions) 
```
dNdxi = [dN1dxi1, dN2dxi1, dN3dxi1;
		 dN1dxi2, dN2dxi2, dN3dxi2]
```
and
```
dNdx = [dN1dx1, dN2dx1, dN3dx1;
        dN1dx2, dN2dx2, dN3dx2]
```
This is "opposite" of the standard matrix if we have ``\partial N_i/\partial x_j``, i.e. the index ``i`` gives the column and ``j`` the row. Therefore, we must calculate the mapping of the ``i``th `MATLAB` as
```
dNdx(:, i) = J' \ dNdxi(:, i)
```
However, to utilize fast vector-matrix multiplication in `MATLAB`, we can simply do all at once,
```
dNdx = J' \ dNdxi
```
which is equivalent to 
```
for i = 1:size(dNdxi, 2)
	dNdx(:, i) = J' \ dNdxi(:, i)
end
```
but faster.

!!! note "Index notation"
    The jacobian ``\underline{\underline{J}} := \partial \underline{x}/\partial \underline{\xi}`` becomes a matrix, whose components are ``J_{ij} = \partial x_i / \partial \xi_j``:
    ```math
		\underline{\underline{J}} = \begin{bmatrix} J_{11} & J_{12} \\ J_{21} & J_{22} \end{bmatrix} = 
		\begin{bmatrix} 
			\partial x_1/\partial \xi_1 & \partial x_1/\partial \xi_2 \\ 
			\partial x_2/\partial \xi_1 & \partial x_2/\partial \xi_2 
		\end{bmatrix}
	```
	If we have a collection of vectors, e.g. ``\underline{x}_\alpha``, and would like to denote the index ``i`` of the vector ``\underline{x}_\alpha``, we enclose the vector in brackets before indexing, i.e. ``[\underline{x}_\alpha]_i``

"""

# ╔═╡ 67c6438f-2547-4580-82fe-a3263203e939
md"""
### Dirichlet boundary conditions
The global shape function, ``N_i`` at the node coordinate, ``\underline{x}_j``, is, as in 1D,
```math
N_i(\underline{x}_j) = \left\lbrace \begin{matrix} 1 & i = j \\ 0 & \text{else} \end{matrix} \right.
```
Check this by considering e.g. node number 4 and manipulate the interactive figure before "Parametric elements". This implies, that when we approximate a function as
```math
T(\underline{x}) \approx T_h(\underline{x}) = \sum_{i = 1}^{N_\mathrm{dofs}} N_i(\underline{x}) a_i
```
we have that ``T_h(\underline{x}_j) = a_j``. So to prescribe the value at a node, we simply prescribe the value of the unknown with that dof-number, just as we did in 1D! When given a mesh, we can request the node numbers on a boundary by using `node_sets`. 
We then build up all the constrained dofs into a vector, e.g. `cdofs`, accompagnied by the known values at the constrained dofs, `ac`, such that we can solve for the unknown dof-values, `af`, by doing,
```
fdofs = setdiff((1:ndofs)', cdofs)
af = K(fdofs, fdofs) \ (f(fdofs) - K(fdofs, cdofs) * ac);
```
"""

# ╔═╡ a4d60a06-d262-4cd1-aad2-f19dd221b4df
md"""
### Neumann boundary conditions
While we saw that Dirichlet BCs were basically equivalent to the 1D case, the Neumann BCs are quite a bit more tricky. Notice that we would like to add the contribution,
```math
f_i^\mathrm{NBC} = -\int_{\Gamma_\mathrm{NBC}} N_i(x) q_\mathrm{n}(x)\ \mathrm{d}\Gamma
```
due to the Neumann BCs on the boundary ``\Gamma_\mathrm{NBC}`` to the load vector ``f_i``. Hence, we need to integrate over the boundary facet: In 2d, this is an edge, and in 3d a face, but we will only consider 2d and edges. 

!!! note "Course curriculum"
	The following **derivations** of the mapping of the integration weights on the element facets in the following description, is outside the scope of the course curriculum. However, **using** the results of the derivation to correctly apply Neumann boundary conditions is part of the curriculum.


To get started, let's remind ourselves about the facets on the reference elements,
"""

# ╔═╡ 8cc38e33-13f2-47a9-8b33-66ce560fbb3b
LocalResource(joinpath(@__DIR__, "refshapes_2d.svg"))

# ╔═╡ cad4715c-72ef-452c-bbf4-7bbd6149ef13
md"""
When we want to apply a Neumann BC, we will get a list of facets, these are described by a matrix where each column contains the global element number and the local facet number. By using the element number, we can get the coordinates of the current element, and by getting the local facet number, we know how to map the coordinates from a reference line [-1, 1] into the local coordinates in the reference element. For the triangle we have with ``s = [\xi^\mathrm{line} + 1]/2``,
```math
\underline{\xi}(\xi^\mathrm{line}) = \left\lbrace \begin{matrix}
[1 - s,\ s]^\mathrm{T} & f = 1 \\
[0,\ s]^\mathrm{T} & f = 2 \\
[s,\ 0]^\mathrm{T} & f = 3 \end{matrix} \right., \quad 

```
and for the quadrilateral, with ``s = \xi^\mathrm{line}``,
```math
\underline{\xi}(\xi^\mathrm{line}) = \left\lbrace \begin{matrix}
[\phantom{+}s,\ -1]^\mathrm{T} & f = 1 \\
[\phantom{+}1,\ \phantom{+}s]^\mathrm{T} & f = 2 \\
[-s,\ \phantom{+}1]^\mathrm{T} & f = 3 \\
[-1,\ -s]^\mathrm{T} & f = 4 \\
\end{matrix} \right.
```

This logic is already provided in
* `triangle_facet_coords`
* `quadrilateral_facet_coords`

For 2D, our integral becomes a line integral along the facet (edge), ``\Gamma^f``,
```math
\int_{\Gamma^f} N_i(\underline{x}) q_\mathrm{n}(\underline{x})\ \mathrm{d}\Gamma
```
We would like to translate this into an integral over a reference line [-1, 1] using ``\xi^\mathrm{line}``. This becomes essentially like integration by substitution, but now we have three "levels" (i.e. two transitions): From reference line to reference element (triangle or quad) coordinates, ``\underline{\xi}(\xi^\mathrm{line})``, and from reference element to physical element, ``\underline{x}(\underline{\xi})``. Using that ``\mathrm{d}\Gamma := |\mathrm{d}\underline{x}|``, we obtain,
```math
\begin{align}
\mathrm{d}\underline{x} &= \frac{\partial\underline{x}}{\partial\underline{\xi}}\cdot\frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}} \mathrm{d}\xi^\mathrm{line} = \underline{\underline{J}} \cdot \frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}\mathrm{d}\xi^\mathrm{line} \\
\int_{\Gamma^f} h(\underline{x})\ \mathrm{d}\Gamma 
&%= \int_{-1}^{1} h(\underline{x}(\underline{\xi}(\xi^\mathrm{line}))) \left\vert\frac{\partial \underline{x}}{\partial \underline{\xi}} \cdot \frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}\right\vert\ \mathrm{d}\xi^\mathrm{line} =
= \int_{-1}^{1} h(\underline{x}(\underline{\xi}(\xi^\mathrm{line}))) \left\vert\underline{\underline{J}} \cdot \frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}\right\vert\ \mathrm{d}\xi^\mathrm{line}
\end{align}
```
Using the mapping above, we arrive at the following derivatives,
```math
\underbrace{\frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}} = \left\lbrace \begin{matrix}
[-0.5,\ -0.5]^\mathrm{T} & f = 1 \\
[\phantom{+}0.0,\ \phantom{+}0.5]^\mathrm{T} & f = 2 \\
[\phantom{+}0.5,\ \phantom{+}0.0]^\mathrm{T} & f = 3 \end{matrix} \right.}_{\text{Reference Triangle}},
\quad\quad
\underbrace{\frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}= \left\lbrace \begin{matrix}
[\phantom{-}1,\ \phantom{-}0]^\mathrm{T} & f = 1 \\
[\phantom{-}0,\ \phantom{-}1]^\mathrm{T} & f = 2 \\
[-1,\ \phantom{-}0]^\mathrm{T} & f = 3 \\
[\phantom{-}0,\ -1]^\mathrm{T} & f = 4 \\
\end{matrix} \right.}_{\text{Reference Quadrilateral}}
```
And we thus have the approximation of the integral,
```math
\int_{\Gamma^f} h(\underline{x})\ \mathrm{d}\Gamma \approx \sum_{q=1}^{N_\mathrm{qp}} h(\underline{x}) \left\vert\underline{\underline{J}}(\underline{\xi}(\xi_q^\mathrm{line})) \cdot \frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}\right\vert\ w_q^\mathrm{line}
```
Actually, 
```math
\underline{\underline{J}}(\underline{\xi}(\xi_q^\mathrm{line})) \cdot \frac{\partial \underline{\xi}}{\partial \xi^\mathrm{line}}
```
gives the direction of the edge scaled with its local length relative the reference line. In `BasicFEM`, the logic for calculating this is provided, except that this vector is rotated to give the direction of the outwards pointing normal, but the magnitude is the same. The functions 
* `triangle_facet_weighted_normal`
* `quadrilateral_facet_weighted_normal`
give this weighted normal, ``\underline{n}_\mathrm{w}^f``, such that the integration simply becomes
"""

# ╔═╡ 510466ef-7197-4e08-b915-1a23eb5fd095
md"""
which is what should be implemented in e.g. `linear_triangle_heat_neumann` and `bilinear_quadrilateral_heat_neumann`. 

!!! note "What you need to know about Neumann BC for heat equation in 2D" 
	* How to use Equation [^facetintegration] to implement neumann boundary conditions
	* Explain the purpose of the factor ``\underline{n}_\mathrm{w}^f`` in this formula
	* Reason what value ``\vert\underline{n}_\mathrm{w}^f\vert`` should take if ``\underline{\underline{J}}`` is constant (Hint: consider ``\int_{\Gamma^f} \mathrm{d}\Gamma``)
	* Apply Neumann boundary conditions that you have implemented (described below)


"""

# ╔═╡ 87c10954-05f0-4866-ae45-df3fcc103c01
md"""
In order to apply Neumann BCs to a mesh, we can do the same as we have done for adding contributions to the force vector from each element, except that we now need to consider each facet (edge) that are part of the boundary where we want to apply the boundary flux. Specifically, we have to loop over each facet in our set
```
for i = 1:size(facet_set, 2)
    element_number = facet_set(1, i);
    facet_number = facet_set(2, i);
    enodes = element_nodes(:, element_number);
    coords = node_coordinates(:, enodes);
    fe = neumann_routine(qn, coords, facet_number, nquadpoints_neumann);
    % Assemble fe into global f as usual 
	% Note that fe has the local numbering of the element
	% I.e for a triangle, length(fe) = 3, so same logic as for a regular element
end
```
"""

# ╔═╡ e840274e-902e-4b28-a8bd-96af78380a20
md"""
### Reaction flux
Similar to the 1D case, after solving the equation systems we know ``f_i`` for all ``i``, and we have
```math
f_i = \int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
For the following mesh, we have a Dirichlet BC on ``\Gamma_\mathrm{right}``, and want to calculate the flux through this boundary, i.e.
```math
q_\mathrm{right} = \int_{\Gamma_\mathrm{right}} q_\mathrm{n}\ \mathrm{d}\Gamma
```
"""

# ╔═╡ c21b70de-f36d-4c91-9a29-8204bb0f8cc2
LocalResource(joinpath(@__DIR__, "reaction_flux_corners.svg"))

# ╔═╡ a471df59-3c75-4cfa-8f19-8612b20c59b0
md"""
First, recall that the shape functions fulfill the following property,
```math
\sum_{i = 1}^{N_\mathrm{s}} N_i(\underline{x}) = 1 \quad \forall \underline{x}\in\Omega
```
Using this property, we then have that
```math
q_\mathrm{right} = \int_{\Gamma_\mathrm{right}} q_\mathrm{n}\ \mathrm{d}\Gamma = \sum_{i = 1}^{N_\mathrm{s}} \int_{\Gamma_\mathrm{right}} N_i(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma
```
Furthermore, we also have that 
```math
N_i(\underline{x}) = 0\ \forall\ \underline{x} \text{ on } \Gamma_{\mathrm{right}} \text{ if } i \notin \mathbb{S}_\mathrm{right} = \lbrace 5,10,15,20,25 \rbrace
```
Hence, 
```math
q_\mathrm{right} = \sum_{i \in \mathbb{S}_\mathrm{right}} \int_{\Gamma_\mathrm{right}} N_i(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma
```
We can now split the calculation of the global vector, ``f_i``, into the different boundaries, ``\Gamma_\mathrm{right}`` and ``\Gamma_\mathrm{other} = \Gamma_\mathrm{top} \cup \Gamma_\mathrm{left} \cup \Gamma_\mathrm{bottom}``, resulting in
```math
\begin{align}
f_i &= \int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{other}} N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
- \int_{\Gamma_\mathrm{right}} N_i(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma \\
\int_{\Gamma_\mathrm{right}} N_i(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma &= 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{other}} N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
- f_i
\end{align}
```
We are only interested in the indices ``i \in \lbrace 5, 10, 15, 20, 25 \rbrace``. Hence, we don't have any contributions from ``\Gamma_\mathrm{left}``. However, the boundaries next to ``\Gamma_\mathrm{right}``: ``\Gamma_\mathrm{top}`` and ``\Gamma_\mathrm{bottom}``, contribute to indices 25 and 5, respectively. So we have
```math
\begin{align}
\int_{\Gamma_\mathrm{right}} N_{25}(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma &= 
\int_\Omega N_{25}(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{top}} N_{25}(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
- f_{25} \\
\int_{\Gamma_\mathrm{right}} N_i(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma &= 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - f_i, \quad i\in\lbrace 10, 15, 20\rbrace \\
\int_{\Gamma_\mathrm{right}} N_5(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma &= 
\int_\Omega N_5(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{bottom}} N_5(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
- f_5
\end{align}
```
So we can easily calculate the contributions from nodes 10, 15, and 20, since we calculate 
```math
f_i^\mathrm{known} = \int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega, \quad i \in \lbrace 10, 15, 20 \rbrace
```
during the regular assembly of ``\underline{\underline{K}}`` and ``\underline{f}``. And if we have **Neumann boundary conditions on the neighboring boundaries**, ``\Gamma_\mathrm{top}`` and ``\Gamma_\mathrm{bottom}``, we also calculate the full
```math
\begin{align}
f_5^\mathrm{known} &= \int_\Omega N_5(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{bottom}} N_5(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma\\
f_{25}^\mathrm{known} &= \int_\Omega N_{25}(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_{\Gamma_\mathrm{top}} N_{25}(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
\end{align}
```
such that the reaction flux becomes
```math
q_\mathrm{right} = \sum_{i \in \mathbb{S}_\mathrm{right}} \left[ f_i^\mathrm{known} - f_i \right], \quad \mathbb{S}_\mathrm{right} = \lbrace 5, 10, 15, 20, 25 \rbrace
```
In `MATLAB`, if we want to calculate the reaction flux on the boundary with node numbers `rdofs` (in this example, `rdofs = 5:5:25`)
```
[K, f_known] = calculate_matrix_and_vector(...) % Add all known contributions, 
										   		% i.e. except unknown qn
a(cdofs) = ac; % Set prescribed temperatures
a(fdofs) = K(fdofs, fdofs) \ (f_known(fdofs) - K(fdofs, cdofs) * ac);
f(fdofs) = f_known(fdofs)
f(cdofs) = K(cdofs, :) * a;
q_right = sum(f_known(rdofs) - f(rdofs));
```
"""

# ╔═╡ 11bd2c69-1744-4789-8d0b-62c23ae7b896
md"""
However, if we have **Dirichlet boundary conditions on a neighboring boundary**, ``\Gamma_\mathrm{top}`` or ``\Gamma_\mathrm{bottom}``, we **cannot calculate the exact contributions** to the flux across each boundary by using this method. In this case, we can make some guess. For example, if we have a Dirichlet BC on ``\Gamma_\mathrm{bottom}``, we can approximate
```math
\begin{align}
\int_{\Gamma_\mathrm{right}} N_5(\underline{x}) q_\mathrm{n}\ \mathrm{d}\Gamma &= \frac{1}{2}\left[\int_\Omega N_5(\underline{x})\ h\ \mathrm{d}\Omega - f_5 \right]\\ 
\int_{\Gamma_\mathrm{bottom}} N_5(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma&= \frac{1}{2}\left[\int_\Omega N_5(\underline{x})\ h\ \mathrm{d}\Omega - f_5 \right]
\end{align}
```
or potentially scale with the edge lengths on either domain. However, note that as we refine the mesh, the contributions to the total reaction force from each element will decrease (proportionally to the element size), such that this error will decrease.

!!! note "Alternative calculation"
	It is also possible to simply evaluate the integral
	```math
		q_\mathrm{right} = \int_{\Gamma_\mathrm{right}} q_\mathrm{n}\ \mathrm{d}\Gamma
	```
	By performing the numerical integration over the boundary (similar code as for applying Neumann boundary conditions), and calculating the heat flux from the temperature field as ``q_\mathrm{n} = \underline{q}\cdot\underline{n}`` with ``\underline{q} = -\underline{\underline{D}}\ \mathrm{grad}(T) \approx -\underline{\underline{D}}\ \mathrm{grad}(N_i) a_i``, given the solution ``a_i``. While this will converge to the correct solution as we refine the mesh, for coarser meshes it can lead to rather large errors in the boundary flux.

"""

# ╔═╡ 55570b97-a717-4a7c-89c4-a2ae941a8e32
md"""
## L9: Transient heat flow
So far, we only considered the stationary energy balance, i.e. when ``\dot{e} = 0``. 
### Strong form 
If we want to consider how the temperature changes with time, we need to consider the change in internal energy with time, and consider the complete strong form,
```math
\dot{e} = h - \mathrm{div}(\underline{q})
```
In the linear theory, we will assume that the internal energy is proportional to the temperature, i.e. 
```math
e(T) = \rho c_\mathrm{p} [T - T_\mathrm{ref}] + e_\mathrm{ref}
```
where ``\rho`` is the material density, ``c_\mathrm{p}`` is the heat capacity (per mass at constant pressure), ``T_\mathrm{ref}`` the temperature at which the internal energy is ``e_\mathrm{ref}``. Inserting this into our strong form, yields,
```math
\rho c_\mathrm{p} \dot{T} = h - \mathrm{div}(\underline{q})
```
"""

# ╔═╡ af3043ce-2863-44d1-8e7e-db5465651c52
md"""
### Weak form
From hereon, we simplify perform the same steps as before, i.e. we get the weak form as
```math
\int_\Omega \delta T \rho c_\mathrm{p} \dot{T}\ \mathrm{d}\Omega = \int_\Omega \delta T h\ \mathrm{d}\Omega - \int_\Omega \delta T \mathrm{div}(\underline{q})\ \mathrm{d}\Omega, \quad \forall\ \delta T(\underline{x})
```
with the arbitrary test function ``\delta T(\underline{x})``. We then use the divergence theorem, exactly as before,
```math
\int_\Omega \delta T \mathrm{div}(\underline{q})\ \mathrm{d}\Omega = \int_\Gamma \delta T\ \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma - \int_\Omega \mathrm{grad}(\delta T)\cdot \underline{q}\ \mathrm{d}\Omega
```
To obtain the final weak form including the Neumann boundary conditions as,
```math
 \int_\Omega \delta T \rho c_\mathrm{p} \dot{T}\ \mathrm{d}\Omega - \int_\Omega \mathrm{grad}(\delta T)\cdot \underline{q}\ \mathrm{d}\Omega = 
\int_\Omega \delta T\ h\ \mathrm{d}\Omega - 
\int_\Gamma \delta T\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
"""

# ╔═╡ ef099aa1-dd65-468e-af83-632ef4e0b535
md"""
### FE form
The difference from before, is that we now are looking for the function ``T(\underline{x}, t)``, i.e. that the temperature is a function of both the coordinate, ``\underline{x}``, and the time, ``t``. We introduce the approximations
```math
\begin{align}
\delta T(\underline{x}) &\approx \sum_{i = 1}^{N_\mathrm{s}} N_i(\underline{x}) c_i = N_i(\underline{x}) c_i \\
T(\underline{x}, t) &\approx \sum_{j = 1}^{N_\mathrm{s}} N_i(\underline{x}) a_i(t) = N_i(\underline{x}) a_i(t)
\end{align}
```
Notice that we have that ``N_i(\underline{x})`` is only a function of the coordinates ``\underline{x}``, and the coefficients are now a function of time, i.e. ``a_i(t)``. Consequently, we have that
```math
\dot{T} \approx N_i(\underline{x}) \dot{a}_i
```
Inserting this into our weak form (using ``N_i`` instead of ``N_i(\underline{x})`` for brevity), we get,
```math
 \int_\Omega N_i c_i \rho c_\mathrm{p} N_j \dot{a}_j\ \mathrm{d}\Omega - \int_\Omega \mathrm{grad}(N_i) c_i)\cdot \underline{q}\ \mathrm{d}\Omega = 
\int_\Omega N_i c_i\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i c_i\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
As before, we can factor out ``c_i``, and get ``c_i r_i = 0`` for all values of ``c_i``, such that our final expression is
```math
 \underbrace{\int_\Omega N_i \rho c_\mathrm{p} N_j\ \mathrm{d}\Omega}_{M_{ij}}\ \dot{a}_j + \underbrace{\int_\Omega \mathrm{grad}(N_i))^\mathrm{T} \underline{\underline{D}}\ \mathrm{grad}(N_j) \ \mathrm{d}\Omega}_{K_{ij}}\ a_j = \underbrace{\int_\Omega N_i\ h\ \mathrm{d}\Omega - \int_\Gamma N_i\ q_\mathrm{n}\ \mathrm{d}\Gamma}_{f_i}
```
where we inserted ``\underline{q} = -\underline{\underline{D}}\ \mathrm{grad}(N_j a_j)``. Using the matrices, we simply have
```math
M_{ij} \dot{a}_j + K_{ij} a_j = f_i, \quad \text{or as matrices,}\quad \underline{\underline{M}}\ \underline{\dot{a}} + \underline{\underline{K}}\ \underline{a} = \underline{f}
```
The new matrix, ``\underline{\underline{M}}``, is called the mass matrix and gives a resistance to change in temperature (similar to how the mass gives resistance to change in velocity in mechanical problems).

In order to solve this problem, we must introduce the concept of *time steps*, i.e. we would like *discretize* the time history as a set of discrete time intervals for which we perform calculations, 
```math
t_i = [t_0, t_1, t_2, \cdots, t_n, t_{n+1}, \cdots, t_N]
```
Let us consider the interval ``[t_n, t_{n+1}]``, where we introduce the approximation of the time derivative of ``a_j`` as
```math
\dot{a}_j \approx \frac{a_j(t_{n+1}) - a_j(t_n)}{t_{n+1} - t_n}
```
For such problems, we need to know the initial conditions, i.e. ``a_j(t_0)``. Then, our task is to find the values at the next time step, i.e. we are looking for ``a_{n+1}`` given ``a_n``. The question then becomes for our equation,
```math
M_{ij} \dot{a}_j + K_{ij} a_j = f_i
```
at what time should we evaluate the term ``K_{ij} a_j``? We can choose to do this at ``t = t_n`` or ``t = t_{n+1}``. We can also evaluate this anywhere inbetween, i.e. at ``t = t_n + [t_{n+1}-t_n]\theta``, where ``\theta \in [0,1]``. Introducing this into our equation yields
```math
M_{ij}\Bigg[ \frac{a_j(t_{n+1}) - a_j(t_n)}{t_{n+1} - t_n}\Bigg] + K_{ij} \Bigg[a_j(t_n) + \theta[a_j(t_{n+1}) - a_j(t_n)]\Bigg] = f_i
```
We then reorganize and multiply with ``\Delta t = t_{n+1} - t_n``, to obtain
```math
\left[M_{ij} + \Delta t \theta K_{ij}\right] a_j(t_{n+1})  = f_i + [M_{ij}- \Delta t [1 - \theta] K_{ij}] a_j(t_n)
```
If we write this in matrix notation we have
```math
\left[\underline{\underline{M}} + \Delta t \theta \underline{\underline{K}}\right] \underline{a}(t_{n+1})  = \underline{f} + \left[\underline{\underline{M}}- \Delta t [1 - \theta]\underline{\underline{K}}\right] \underline{a}(t_n)
```
And we can solve for the unknowns at the next time step as
```math
\underline{a}(t_{n+1}) = \left[\underline{\underline{M}} + \Delta t \theta \underline{\underline{K}}\right]^{-1}\left[ \underline{f} + \left[\underline{\underline{M}}- \Delta t [1 - \theta]\underline{\underline{K}}\right] \underline{a}(t_n)\right]
```
Choosing ``\theta = 0`` results in a so-called fully explicit time integration (Forward Euler), which can be faster (since the matrix to be factorized is constant). However, this method requires small time steps to be stable, in general 
```math
\Delta t \leq \frac{2}{[1 - 2\theta]\lambda_\mathrm{max}}, \quad \lambda_\mathrm{max} = \text{ maximum eigenvalue of } \left[\underline{\underline K} - \lambda \underline{\underline M}\right] \underline{\Lambda} = 0
```
is a stable time step. Choosing ``\theta \geq 0.5`` is unconditionally stable. Often, the fully implicit ``\theta = 1`` is chosen. This is called Backward Euler.
"""

# ╔═╡ 509c3307-b814-4da9-aecc-7adb47f8ec95
md"""
## L10: Introduction to Abaqus
See notes on Canvas
## L11: Guest lectures
See presentations uploaded on Canvas
"""

# ╔═╡ 1f5672bc-28af-4ab7-b159-197ebf1e12a3
md"""
### Cauchy's theorem

#### The traction vector
To describe a distributed load acting on a surface, we need to define the traction, ``\underline{t}``, as the force per area. In general, the traction is not constant, and is defined as the ratio between the force, ``\underline{F}``, acting on a small patch with area ``A``, as we let the patch size go to zero, specifically
```math
\underline{t} := \lim_{A\rightarrow 0} \frac{\underline{F}}{A}
```
The surface that we consider can be either a physical surface of a body, or we can create an imaginary surface inside a continuum by cutting it and considering the internal forces acting on this surface. The latter is illustrated by cutting the potato-shaped body below with a plane having the normal vector $\underline{n}$.
To describe a distributed load acting on a surface, we need to define the traction, ``\underline{t}``, as the force per area. In general, the traction is not constant, and is defined as the ratio between the force, ``\underline{F}``, acting on a small patch with area ``A``, as we let the patch size go to zero, specifically
```math
\underline{t} := \lim_{A\rightarrow 0} \frac{\underline{F}}{A}
```
The surface that we consider can be either a physical surface of a body, or we can create an imaginary surface inside a continuum by cutting it and considering the internal forces acting on this surface. The latter is illustrated by cutting the potato-shaped body below with a plane having the normal vector ``\underline{n}``.
"""

# ╔═╡ 4167e37d-b9e2-484f-935c-729e0507630b
LocalResource(joinpath(@__DIR__, "traction_definition.svg"))

# ╔═╡ bebe6b31-f5f3-4719-b970-7f0583bc3674
md"""
#### Multiple traction vectors at the same point
When using the traction on a specific plane to investigate the load inside a body, as shown above, we only evaluate the load on a specific plane. But at the same point in the body, we could cut it using different planes, as illustrated below with the normal vectors ``\underline{n}_1`` and ``\underline{n}_2``.
"""

# ╔═╡ 0d31dfd6-df4b-4d08-a099-c69c2f2cb071
LocalResource(joinpath(@__DIR__, "traction_multiple_cuts.svg"))

# ╔═╡ 9459013d-e12c-4a82-9636-5ce912cccaf3
md"""
This leads to different traction vectors, ``\underline{t}_1`` and ``\underline{t}_2``, which is inconvenient if we want to check if the material will fail in this location. Then we would need to calculate and check the traction on all planes.

In this lecture, we aim to find a quantity that can describe the loading on all planes in the body, instead of having a traction vector for each plane. To do this, let's cut out a small tetrahedron at the point of interest.
"""

# ╔═╡ 48ebeebf-8507-4941-9361-4a5b4fe601db
begin
	σtet_L1_slider = @bind σtet_L1 Slider(0.1:0.1:1.0; default=0.5)
	σtet_L2_slider = @bind σtet_L2 Slider(0.1:0.1:1.0; default=0.5)
	σtet_L3_slider = @bind σtet_L3 Slider(0.1:0.1:1.0; default=0.5)
	σtet_azimuth_slider = @bind σtet_azimuth Slider(range(-π, π, 24+1); default=-3π/12)
	σtet_elevation_slider = @bind σtet_elevation Slider(range(-π, π, 24+1); default=π/6)
	md"""
	#### The stress tetrahedron
	We consider a tetrahedron, where 3 sides are aligned with the coordinate system as
	
	| Color | Plane | Normal | Area |
	| ----- | ----- | ------ | ---- |
	| $\textcolor{blue}{\text{Blue}}$  | $X_2$-$X_3$ | $\boldsymbol{n}_1=-\boldsymbol{e}_1$ | $A_1$ |
	| $\textcolor{orange}{\text{Orange}}$   | $X_1$-$X_3$ | $\boldsymbol{n}_2=-\boldsymbol{e}_2$ | $A_2$ |
	| $\textcolor{green}{\text{Green}}$ | $X_1$-$X_2$ | $\boldsymbol{n}_3=-\boldsymbol{e}_3$ | $A_3$ |
	| $\textcolor{cyan}{\text{Cyan}}$ |   | $\hat{\boldsymbol{n}}$ | $\hat{A}$ |
	
	Defining the lengths, $L_1$, $L_2$, and $L_3$, along the coordinate axes fully defines the tetrahedron. 
	
	The traction vectors acting on each side, $\boldsymbol{t}$, are shown in grey. You can change the side lengths, and rotate the figure with the following controls to better understand the geometry.
	
	|  |  |  |
	|---|---|---|
	|``L_1``: $(σtet_L1_slider) | ``L_2``: $(σtet_L2_slider) | ``L_3``: $(σtet_L3_slider) |
	| Rotation: $(σtet_azimuth_slider) | Elevation: $(σtet_elevation_slider) | | 
	"""
end

# ╔═╡ a9e3d9df-5ab5-4ced-a4e6-c61ce2e79046
begin
	# utils
	mkobsvec(n) = Plt.Observable(fill(NaN, n))
	convert_to_point(v::Vec{dim}) where dim = convert(Plt.Point{dim}, v)
	convert_to_vec(v::Vec{dim}) where dim = convert(Plt.Vec{dim}, v)
	convert_to_point(v::NTuple{dim}) where dim = convert(Plt.Point{dim}, v)
	function calculate_relative_areas(L1, L2, L3)
		n0 = zero(Vec{3})
		n1 = Vec{3}(( L1, 0.0, 0.0))
		n2 = Vec{3}((0.0,  L2, 0.0))
		n3 = Vec{3}((0.0, 0.0,  L3))

		faces = [(n0, n3, n2), (n0, n1, n3), (n0, n2, n1), (n1, n2, n3)]
		face_areas = map(face -> norm((face[2] - face[1]) × (face[3] - face[1])), faces)
		return face_areas[1:3] / face_areas[4]
	end
	function calculate_quantities_stress_tetrahedron(lengths = (0.8, 1.2, 1.5))
		L1, L2, L3 = lengths
		n_scaleface = 0.25
		n0 = zero(Vec{3})
		n1 = Vec{3}(( L1, 0.0, 0.0))
		n2 = Vec{3}((0.0,  L2, 0.0))
		n3 = Vec{3}((0.0, 0.0,  L3))
		
		σ = SymmetricTensor{2,3}((1.0, 0.5, 0.0, 0.5, 0.2, 0.2))/2
		
		faces = [(n0, n3, n2), (n0, n1, n3), (n0, n2, n1), (n1, n2, n3)]
		face_colors = [(:blue, 1.0), (:orange, 1.0), (:green, 1.0), (:cyan, 0.3)]
		face_colors = map(first, face_colors)
		v_tmp = Float64[]
		v_colors = []
		face_normals = []
		face_tractions = []
		face_centers = []
		face_areas = Float64[]
		for (face, color) in zip(faces, face_colors)
			push!(face_centers, sum(face)/length(face))
			w = (face[2] - face[1]) × (face[3] - face[1])
			n = w / norm(w)
			push!(face_tractions, n ⋅ σ)
			push!(face_areas, norm(w))
			push!(face_normals, n_scaleface * n)
			for n in face
				push!(v_colors, color)
				append!(v_tmp, n...)
			end
		end
		vertices= collect(transpose(reshape(v_tmp, 3, :)))
		face_vertices = collect(transpose(reshape(1:12, 3, :)))
		return vertices, face_vertices, v_colors, face_colors, map(convert_to_point, face_centers), map(convert_to_vec, face_normals), map(convert_to_vec, face_tractions)
	end

	function update_stress_tetrahedron!(obs, L1, L2, L3, azimuth, elevation)
		verts, f_verts, v_colors, f_colors, f_centers, f_normals, f_tractions = calculate_quantities_stress_tetrahedron((L1, L2, L3))
		obs["verts"][] = verts
		obs["f_centers"][] = f_centers
		obs["f_normals"][] = f_normals
		obs["f_tractions"][] = f_tractions
		obs["azimuth"][] = azimuth
		obs["elevation"][] = elevation
		return nothing
	end
	
	function setup_stress_tetrahedron()
		fig = Plt.Figure(size=(600,400))
		azimuth = Plt.Observable(-0.25π)
		elevation = Plt.Observable(π/8)
		ax = Plt.Axis3(fig[1,0:2]; xlabel=L"X_1", ylabel=L"X_2", zlabel=L"X_3", azimuth, elevation)
		foreach(lf -> lf(ax, (-0.3, 1.1)), (Plt.xlims!, Plt.ylims!, Plt.zlims!))
		
		# Strategy is to plot each face with separate nodes, so 3*4 nodes in total
		verts, f_verts, v_colors, f_colors, f_centers, f_normals, f_tractions = calculate_quantities_stress_tetrahedron((0.8, 1.2, 1.5))
		obs = Dict(
			"verts" => Plt.Observable(verts),
			"f_centers" => Plt.Observable(f_centers),
			"f_normals" => Plt.Observable(f_normals),
			"f_tractions" => Plt.Observable(f_tractions),
			"azimuth" => azimuth, "elevation" => elevation
		)
		
		Plt.mesh!(ax, obs["verts"], f_verts; color=v_colors)
		
		p_normals = Plt.arrows2d!(ax, 
								  obs["f_centers"], 
								  obs["f_normals"]; 
								  color=f_colors,
								  #linewidth=0.035, 
								  #arrowsize=Plt.Vec3f(0.10, 0.10, 0.15)
								 )
		p_traction = Plt.arrows2d!(ax, 
								   obs["f_centers"], 
								   obs["f_tractions"]; 
								   color=(:grey, 1.0),
								   #linewidth=0.040, 
								   #arrowsize=Plt.Vec3f(0.10, 0.10, 0.15)
								   )

		zeropoints = [zero(Plt.Point3) for _ in 1:3]
		endpoints = [Plt.Point3(ntuple(i->1.0*(i==j), 3)) for j in 1:3]
		Plt.hidespines!(ax)
		Plt.hidedecorations!(ax)
		
		Plt.arrows2d!(ax, 
					  zeropoints, 
					  endpoints; 
					  #linewidth=0.02, 
					  #arrowsize=Plt.Vec3f(0.075, 0.075, 0.1)
					  )
		Plt.text!(ax, [Plt.Point3(1.0, 0, 0.05), Plt.Point3(0, 1.0, 0.05), Plt.Point3(0.00, 0.00, 1.10)]; text=[L"\mathbf{e}_1", L"\mathbf{e}_2", L"\mathbf{e}_3"], fontsize=20)

		arrowpoints = let # local scope
			w_tail = 0.2; l_tail = 0.6; w_head = 0.5
			l_head = 1 - l_tail
			Plt.Point2f[(0, -w_tail/2), (l_tail, -w_tail/2), (l_tail, -w_head/2), (l_tail + l_head, 0.0), (l_tail, w_head/2), (l_tail, w_tail/2), (0, w_tail/2)] .+ Plt.Point2f(0, 0.5)
		end

		face_legend_elements = [Plt.PolyElement(;color, points=arrowpoints) for color in f_colors]
		face_legends = [L"\mathbf{n}_1", L"\mathbf{n}_2", L"\mathbf{n}_3", L"\hat{\mathbf{n}}"]
		traction_legend_element = Plt.PolyElement(;color=:grey, points=arrowpoints)
		
		Plt.Legend(fig[0,1], [face_legend_elements..., traction_legend_element], [face_legends..., L"\mathbf{t}"];
		patchsize = (20, 15), orientation = :horizontal, labelsize=20,
		tellwidth=true, tellheight=true)
		
		return fig, obs
	end
	stress_tetrahedron_fig, stress_tetrahedron_obs = setup_stress_tetrahedron();
end;

# ╔═╡ 75a525f7-e727-4c55-9fcd-a4335d97c432
let
	update_stress_tetrahedron!(stress_tetrahedron_obs, σtet_L1, σtet_L2, σtet_L3, σtet_azimuth, σtet_elevation)
	stress_tetrahedron_fig
end

# ╔═╡ dc02ad0d-058b-4369-be78-44039aa6e6b1
md"""
#### Force equilibrium 
We require that the tetrahedron is in equilibrium, i.e.
```math
\sum \underline{F} = \underline{t}_1 A_1 + \underline{t}_2 A_2 + \underline{t}_3 A_3 + \hat{\underline{t}} \hat{A} + \underline{b} V = \underline{t}_i A_i + \hat{\underline{t}} \hat{A}  + \underline{b} V = 0
```
where ``\underline{b}`` is a body load (force per volume) and ``V=A_3 L_3 / 3`` the volume of the tetrahedron. Dividing by ``A_3`` and letting the size go towards zero, we see that the body load term vanishes in comparison to the other terms, and in the limit as the size goes to zero, we only require that 
```math
\underline{t}_i A_i + \hat{\underline{t}} \hat{A} = 0
```
To move forward, we need to find a relationship between the areas ``A_i`` and ``\hat{A}``
"""

# ╔═╡ 3345a144-6871-4b5b-b60b-db702770b30a
md"""
#### Area relationships
To find a relationship between the areas, we will use Gauss' divergence theorem,
```math
\int_\Omega \mathrm{div}(\underline{a})\ \mathrm{d}\Omega = \int_\Gamma \underline{a}\cdot\underline{n}\ \mathrm{d}\Gamma
```
where ``\underline{n}`` is the outwards-pointing normal vector from the body ``\Omega`` with boundary ``\Gamma``. ``\underline{a}`` a vector-field defined in ``\Omega``. To study the areas, we will choose the spatially constant base vector ``\underline{e}_i`` as our vector field, yielding
```math
0 = \int_\Gamma \underline{e}_i \cdot \underline{n}\ \mathrm{d}\Gamma = 
\int_\Gamma n_i\ \mathrm{d}\Gamma, \quad \Rightarrow \quad \int_\Gamma \underline{n}\ \mathrm{d}\Gamma = \underline{0}
```
as ``\mathrm{div}(\underline{e}_i) = 0``. 
As seen above, our tetrahedron consists of 4 flat sides, allowing us to write the integral as a sum over all sides, i.e. 
```math
\int_\Gamma \underline{n}\ \mathrm{d}\Gamma = \hat{A}\ \hat{\underline{n}} + A_i \underline{n}_i = \underline{0}
```
Using that ``\underline{n}_i = -\underline{e}_i``, we get ``\hat{A}\ \hat{\underline{n}} = A_i \underline{e}_i``. If we take the dot product with the base vector ``\underline{e}_j``, we get
```math
\begin{align}
\hat{A}\ \hat{\underline{n}} \cdot \underline{e}_j &= A_i \underline{e}_i \cdot \underline{e}_j = A_i \delta_{ij} \\
\hat{A}\ \hat{n}_j &= A_j
\end{align}
```

"""

# ╔═╡ b5815863-6379-4711-a9fd-4cc69e8388e1
let
	A1, A2, A3 = calculate_relative_areas(σtet_L1, σtet_L2, σtet_L3)
	fmt = FormatSpec("0.2f")
	md"""
	For the current geometry above, we have the normal vector ``\hat{\underline{n}} = `` ($(pyfmt(fmt, A1)), $(pyfmt(fmt, A2)), $(pyfmt(fmt, A3))), and the area relations become
	
	``A_1`` = $(pyfmt(fmt, A1)) ``\hat{A}``, 
	``\quad`` 
	``A_2`` = $(pyfmt(fmt, A2)) ``\hat{A}``,
	``\quad`` 
	``A_3`` = $(pyfmt(fmt, A3)) ``\hat{A}``

	We have thus managed to express the areas $A_i$ relative $\hat{A}$ from the normal vector $\hat{\underline{n}}$.
	"""
end

# ╔═╡ 2ec78a16-f64e-4795-8305-d1a36899eb0a
md"""
#### The Cauchy stress
To introduce the Cauchy stress, let's evaluate the equilibrium equation above in each coordinate direction, by taking the dot product with the base vector ``\underline{e}_j``,
```math
A_i \underline{t}_i \cdot \underline{e}_j   + \hat{A} \hat{\underline{t}} \cdot \underline{e}_j = 0
```
Next, we insert the area relationship, ``A_i = \hat{A}\ \hat{n}_i``, and divide by ``\hat{A}``, to get
```math
\hat{n}_i \underline{t}_i \cdot \underline{e}_j + \hat{\underline{t}} \cdot \underline{e}_j = 0
```
Last, we define the Cauchy stress, ``\sigma_{ij} := - \underline{t}_i \cdot \underline{e}_j``, resulting in 
```math
\hat{t}_j = \hat{n}_i \sigma_{ij}
```
In matrix-vector notation, this becomes
```math
\hat{\underline{t}} = \hat{\underline{n}} \cdot \underline{\underline{\sigma}}
```
Since the plane with normal ``\hat{\underline{n}}`` is an arbitrary plan, this result, called the Cauchy stress theorem, shows that the quantity, ``\underline{\underline{\sigma}}``, which is the Cauchy stress, allows us to describe the traction ``\hat{\underline{t}}`` on any plane. I.e., the stress matrix ``\underline{\underline{\sigma}}``, fully describes the load on the material at a point in the body. We only used the ``\hat{\underline{n}}`` and ``\hat{\underline{t}}`` in this derivation, and we will state the Cauchy's theorem as
```math
t_j = n_i \sigma_{ij}, \quad \text{ or in matrix-vector}, \quad \underline{t} = \underline{n} \cdot \underline{\underline{\sigma}} = \underline{n}^\mathrm{T}\underline{\underline{\sigma}}
```
where ``\underline{t}`` is the traction vector on a plane with normal vector ``\underline{n}``. 
"""

# ╔═╡ 6ccce575-56b5-4fc6-969e-b71348c44a81
md"""
### Translational equilibrium

Starting from Newton's 2nd law for a particle, ``\underline{F} = m\underline{a} = m\underline{\ddot{u}}`` (or ``F_j = m \ddot{u}_j``), we can set up this law for a body (or part of a body), ``\Omega`` with boundary ``\Gamma`` as 
```math
\underbrace{\int_\Gamma t_j\ \mathrm{d}\Gamma + \int_\Omega b_j\ \mathrm{d}\Omega}_{F_j} = \underbrace{\int_\Omega \rho \ddot{u}_j\ \mathrm{d}\Omega}_{m\ddot{u}_j}
```
where ``t_j`` is the traction vector on the boundary ``\Gamma``, ``b_j`` is the body (volume) force, ``\rho`` the density and ``\ddot{u}_j`` the acceleration (2nd time-derivative of the displacements ``u_j``). Inserting Cauchy's theorem, we get
```math
\int_\Gamma n_i \sigma_{ij}\ \mathrm{d}\Gamma + \int_\Omega b_j\ \mathrm{d}\Omega = \int_\Omega \rho \ddot{u}_j\ \mathrm{d}\Omega
```
Then we can apply the divergence theorem ,
```math
\int_\Gamma n_i \sigma_{i\textcolor{red}{j}}\ \mathrm{d}\Gamma = \int_\Omega \frac{\partial \sigma_{i\textcolor{red}{j}}}{\partial x_i}\ \mathrm{d}\Omega
```
To see that this is the same as for the case of a vector field, Equation [^divergencetheorem], notice that this is the same as that equation for each distinct value of the the index ``\textcolor{red}{j}``. Including this gives,
```math
\int_\Omega \frac{\partial \sigma_{ij}}{\partial x_i}\ \mathrm{d}\Omega + \int_\Omega b_j\ \mathrm{d}\Omega = \int_\Omega \rho \ddot{u}_j\ \mathrm{d}\Omega
```
Which allows us to apply the localization argument (must hold for any subdomain ``\Omega``, hence it needs to hold point-wise), i.e. 
```math
\frac{\partial \sigma_{ij}}{\partial x_i} + b_j = \rho \ddot{u}_j
```
"""

# ╔═╡ 20f8718e-60ea-415b-aab3-17a90c7f31c6
md"""
### Rotational equilibrium
We consider a box with side lengths ``\Delta x_1 \rightarrow 0`` and ``\Delta x_2 \rightarrow 0``, with a stress ``\underline{\underline{\sigma}}``, and take the torque balance around the center of the box:
"""

# ╔═╡ d7e7d8c5-0673-4c78-a484-4ee8116f1a76
LocalResource(joinpath(@__DIR__, "stress_symmetry_viapdf.svg"))

# ╔═╡ 7cec58dd-dd80-4ce5-b885-f1e2527e481d
md"""
This torque balance around the center of the box, with positive direction counter-clockwise, considering a thickness ``\Delta x_3``, becomes
```math
\begin{align}
\sum T = 0 &= 
\frac{\Delta x_1}{2} \left[\Delta x_2 \Delta x_3 \left[\underline{t}_1^+\right]_2\right] - 
\frac{\Delta x_2}{2} \left[\Delta x_1 \Delta x_3 \left[\underline{t}_2^+\right]_1\right]\\ 
&- 
\frac{\Delta x_1}{2} \left[\Delta x_2 \Delta x_3 \left[\underline{t}_1^-\right]_2\right] + 
\frac{\Delta x_2}{2} \left[\Delta x_1 \Delta x_3 \left[\underline{t}_2^-\right]_1\right]
\end{align}
```
Dividing by ``\Delta x_1 \Delta x_2 \Delta x_3 / 2`` and inserting the components of ``\underline{t}_i^\pm``, gives
```math
\begin{align}
0 &= 
\left[\underline{t}_1^+\right]_2 - 
\left[\underline{t}_2^+\right]_1 - 
\left[\underline{t}_1^-\right]_2 + 
\left[\underline{t}_2^-\right]_1 \\
0 &= \sigma_{12} - \sigma_{21} - [-\sigma_{12}] + [-\sigma_{21}] = 2\sigma_{12} - 2\sigma_{21} \Rightarrow \sigma_{12} = \sigma_{21}
\end{align}
```
Showing that the stress (in 2D) is symmetric. The same results follows by doing the same derivation in 3D. 
"""

# ╔═╡ f13237a3-1bf4-4ec9-a797-20458df04e52
md"""
## L12b: Linear elasticity
Hooke's law using the stiffness tensor, ``\underline{\underline{D}}``, in Voigt notation
```math
\underline{\sigma} = \underline{\underline{D}}\ \underline{\varepsilon}
```
with the strain being
```math
\underbrace{\underline{\varepsilon}}_{3\text{D}} = \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\ 2\varepsilon_{23} \\ 2\varepsilon_{13} \\ 2\varepsilon_{12} \end{bmatrix}
, \quad
\underbrace{\underline{\varepsilon}}_{2\text{D}} = \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{bmatrix}
```
with the components ``\varepsilon_{ij}`` defined as
```math
\varepsilon_{ij} = \frac{1}{2}\left[\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right]
```

### 3D case
```math
\underline{\underline{D}} = \frac{E}{[1+\nu][1 - 2\nu]} \begin{bmatrix}
1 - \nu & \nu & \nu & 0 & 0 & 0 \\ 
\nu & 1 - \nu & \nu & 0 & 0 & 0 \\
\nu & \nu & 1 - \nu & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{1 - 2\nu}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \frac{1 - 2\nu}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \frac{1 - 2\nu}{2}
\end{bmatrix}
```
### Plane strain
"Plane strain" implies that we assume that ``\varepsilon_{ij} = 0`` if ``i = 3`` or ``j = 3``. In this case, the stiffness tensor in Voigt format becomes
```math
\underline{\underline{D}} = \frac{E}{[1+\nu][1 - 2\nu]} \begin{bmatrix}
1 - \nu & \nu  & 0 \\ 
\nu & 1 - \nu  & 0 \\
0 & 0 & \frac{1 - 2\nu}{2}
\end{bmatrix}
```

The plane strain assumption holds if the out-of-plane motion is constrained, and the loading is uniform across the thickness. Note that this is always an approximation! The typical example is when the thickness (out-of-plane dimension) is large, but this doesn't always hold. Hence the key question to ask is if the out-of-plane motion of the analyzed cross-section is contrained.  

!!! note "Out-of-plane stress components"
    When calculating the stress under plane strain, we get a 2d stress (3 components
	in the Voigt format), which is enough for ensuring equilibrium. However, since ``\sigma_{ij} \neq 0`` for ``i = 3`` or ``j = 3``, we have to include those stress components when calculating effective stress measures. An easy way to obtain the full stress matrix, is using the stiffness tensor for 3d, and inserting zeros for the strain components that should be zero. 

### Plane stress
"Plane stress" implies that we assume that ``\sigma_{ij} = 0`` if ``i = 3`` or ``j = 3``. In this case, the stiffness tensor in Voigt format becomes
```math
\underline{\underline{D}} = \frac{E}{1-\nu^2} \begin{bmatrix}
1 & \nu  & 0 \\ 
\nu & 1  & 0 \\
0 & 0 & \frac{1 - \nu}{2}
\end{bmatrix}
```

The plane stress assumption holds if the out-of-plane motion is unconstrained, and the loading is uniform across the thickness. The typical example is when the thickness (out-of-plane dimension) is small, but there are cases when plane stress is a good assumption also for a large of out plane dimensions. Hence the key question to ask is if the out-of-plane motion of the analyzed cross-section is unconstrained. Note that for thin cross-sections, we may get a good approximation even if the loading is not uniform across the thickness direction. 

!!! note "Out-of-plane strain components"
    When calculating the stress under plane strain, we use 2d strain (3 components
	in the Voigt format), which is enough for ensuring equilibrium. If we need to obtain the nonzero strain components for ``i = 3`` or ``j = 3``, we can use the the stiffness tensor for 3d, and inserting zeros for the stress components that should be zero, before calculating ``\underline{\varepsilon} = \underline{\underline{D}}^{-1} \underline{\sigma}``. 

"""

# ╔═╡ e2f98607-c9a6-4a89-ba81-87b5a83cb789
md"""
### Effective stress measures and yielding
Using *linear elasticity* is an assumption, which gives a very good approximation for many engineering materials when the stress is not too high. The question of what is a too high stress is often answered by calculating effective stress measures. For metals, we often want to stay below the yield limit of the material (after which we induce plastic deformations), and a typical effective stress measure is the **von Mises** effective stress, ``\sigma_\mathrm{vM}``,
```math
\sigma_\mathrm{vM} = \sqrt{\frac{
[\sigma_{11} - \sigma_{22}]^2 + [\sigma_{11} - \sigma_{33}]^2 + [\sigma_{22} - \sigma_{33}]^2 + 6[\sigma_{12}^2 + \sigma_{13}^2 + \sigma_{23}^2]}{2}}
```
For metals, the value of ``\sigma_\mathrm{vM}`` can be compared to the yield limit, to see if it is valid to consider linear elasticity and if no permanent deformations will occur. What is special about the von Mises effective stress is that it does not depend on the so-called hydrostatic stress, i.e. ``\sigma_\mathrm{vM}(\underline{\underline{\sigma}}) = \sigma_\mathrm{vM}(\underline{\underline{\sigma}} - p \underline{\underline{I}})``, where ``p`` is the pressure. This is a good assumption for metals, but for granular materials, such as soil, gravel, and concrete, the **Drucker-Prager** criterion may be more suitable. It is defined as
```math
\sigma_\mathrm{DP} = \sigma_\mathrm{vM} - B\underbrace{[\sigma_{11} + \sigma_{22} + \sigma_{33}]}_{-p}
```
These are just some examples of different stress measures, and what to use really depends on the specific analysis to be made. While these are suitable for e.g. determining the maximum load on a structure, they are **not suitable** for analyzing fatigue. 
"""

# ╔═╡ 5f0cf5a1-bd8c-4637-a28d-e4fd134254ab
md"""
## L13: Weak form of the mechanical equilibrium
```math
\int_\Omega \delta u_j\ \frac{\partial \sigma_{ij}}{\partial x_i}\ \mathrm{d}\Omega + \int_\Omega \delta u_j\ b_j\ \mathrm{d}\Omega = \int_\Omega \delta u_j\ \rho \ddot{u}_j\ \mathrm{d}\Omega
```
Differentiate ``\partial[\delta u_j \sigma_{ij}] / \partial x_i``,
```math
\frac{\partial}{\partial x_i}[\delta u_j \sigma_{ij}] = \frac{\partial \delta u_j}{\partial x_i} \sigma_{ij} + \delta u_j \frac{\partial \sigma_{ij}}{\partial x_i}
```
Denoting ``v_i = \delta u_j \sigma_{ij}``, we can apply the divergence theorem for a vector, i.e. 
```math
\int_\Omega \frac{\partial}{\partial x_i}[\delta u_j \sigma_{ij}]\ \mathrm{d}\Omega = \int_\Omega \frac{\partial v_i}{\partial x_i}\ \mathrm{d}\Omega = \int_\Gamma n_i v_i\ \mathrm{d}\Gamma = \int_\Gamma \delta u_j n_i \sigma_{ij}\ \mathrm{d}\Gamma = \int_\Gamma \delta u_j t_j\ \mathrm{d}\Gamma
```
And putting it all together, we get
```math
\int_\Gamma \delta u_j t_j\ \mathrm{d}\Gamma - \int_\Omega \frac{\partial \delta u_j}{\partial x_i} \sigma_{ij}\ \mathrm{d}\Omega + \int_\Omega \delta u_j\ b_j\ \mathrm{d}\Omega = \int_\Omega \delta u_j\ \rho \ddot{u}_j\ \mathrm{d}\Omega
```
And reorganizing we get,
```math
\int_\Omega \delta u_j\ \rho \ddot{u}_j\ \mathrm{d}\Omega + \int_\Omega \frac{\partial \delta u_j}{\partial x_i} \sigma_{ij}\ \mathrm{d}\Omega = \int_\Gamma \delta u_j t_j\ \mathrm{d}\Gamma + \int_\Omega \delta u_j\ b_j\ \mathrm{d}\Omega
```
We note that the stress, ``\sigma_{ij}``, and the gradient ``\partial \delta u_j/\partial x_i``, are 2nd-order objects (i.e. two indices). We can represent these as matrices, e.g.
```math
\underline{\underline{\sigma}} = \begin{bmatrix} 
\sigma_{11} & \sigma_{12} \\ \sigma_{21} & \sigma_{22} \end{bmatrix}
```
And we define a special multiplication using ``:``, as 
```math
\begin{align}
\underline{\underline{a}}:\underline{\underline{b}} &= a_{ij} b_{ij}\\
a_{ij} b_{ij} &= a_{11} b_{11} + a_{21} b_{21} + a_{12} b_{12} + a_{22} b_{22}, \quad \text{in 2d}
\end{align}
```
We can then write the weak form as, using that ``\sigma_{ij} = \sigma_{ji}``,
```math
\int_\Omega \delta \underline{u}\ \cdot \rho \ddot{\underline{u}}\ \mathrm{d}\Omega + \int_\Omega \frac{\partial \delta \underline{u}}{\partial \underline{x}} :\underline{\underline{\sigma}}\ \mathrm{d}\Omega = \int_\Gamma \delta \underline{u} \cdot \underline{t}\ \mathrm{d}\Gamma + \int_\Omega \delta \underline{u}\cdot\underline{b}\ \mathrm{d}\Omega
```
"""

# ╔═╡ 976fbb59-1adf-4e95-84b8-75c7bf302917
md"""
## L14: Finite element analysis of linear elasticity
We now introduce the approximations,
```math
\begin{align}
\delta\underline{u}(\underline{x}) &\approx \sum_{i = 1}^{N_\mathrm{s}} \underline{M}_i(\underline{x}) c_i = \underline{M}_i(\underline{x}) c_i \\
\underline{u}(\underline{x}) &\approx \sum_{j = 1}^{N_\mathrm{s}} \underline{M}_j(\underline{x}) a_j = \underline{M}_j(\underline{x}) a_j
\end{align}
```
where ``\underline{M}_i(\underline{x})`` are *vector-valued* shape functions. We will get back to how exactly these are defined later

!!! note "Notation"
    The same symbol, ``N_i(\underline{x})`` and ``\underline{N}_i(\underline{x})``, 
    is used for both scalar and vector-valued shape functions. Here we use
    ``N_i(\underline{x})`` and ``\underline{M}_i(\underline{x})`` to 
    clearly separate the two).

Inserting this into the weak form, we obtain (skipping ``(\underline{x})`` in ``\underline{M}_i(\underline{x})`` for brevity),
```math
\int_\Omega [\underline{M}_i c_i]\cdot \rho \underline{M}_j \ddot{a}_j\ \mathrm{d}\Omega + \int_\Omega \frac{\partial [\underline{M}_i c_i]}{\partial \underline{x}} :\underline{\underline{\sigma}}\ \mathrm{d}\Omega = \int_\Gamma [\underline{M}_i c_i] \cdot \underline{t}\ \mathrm{d}\Gamma + \int_\Omega [\underline{M}_i c_i]\cdot\underline{b}\ \mathrm{d}\Omega
```
And just as for the 1D and 2D heat equation, we can factor out the arbitrary coefficients, ``c_i``, to obtain ``c_i r_i = 0``, resulting in that we require ``r_i = 0``, and we get 
```math
\int_\Omega \underline{M}_i\cdot \rho \underline{M}_j\ \mathrm{d}\Omega\ \ddot{a}_j + \int_\Omega \frac{\partial \underline{M}_i}{\partial \underline{x}} :\underline{\underline{\sigma}}\ \mathrm{d}\Omega = \int_\Gamma \underline{M}_i \cdot \underline{t}\ \mathrm{d}\Gamma + \int_\Omega \underline{M}_i\cdot\underline{b}\ \mathrm{d}\Omega
```
While it is possible to work with this form, many Finite Element codes use the *Voigt* format to represent the matrices, ``\partial \underline{M}_i/\partial \underline{x}`` and ``\underline{\underline{\sigma}}``, as vectors. Specifically we define the Voigt-representation of the stress, ``\underline{\sigma}``, as
```math
\underbrace{\underline{\sigma}}_{3\text{D}} = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix}
, \quad
\underbrace{\underline{\sigma}}_{2\text{D}} = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix}
```
where we only store 6 or 3 components in 3d and 2d, instead of the full 9 or 4 components due to the symmetry, ``\sigma_{ij} = \sigma_{ji}``. 
The gradient of the shape functions, ``\partial \underline{M}_i/\partial \underline{x}``, is represented by the B-vector, ``\underline{B}_i``,

```math
\underbrace{\underline{B}_i}_{3\text{D}} = \begin{bmatrix} 
\frac{\partial [\underline{M}_i]_1}{\partial x_1} \\ 
\frac{\partial [\underline{M}_i]_2}{\partial x_2} \\ 
\frac{\partial [\underline{M}_i]_3}{\partial x_3} \\ 
\frac{\partial [\underline{M}_i]_2}{\partial x_3} + \frac{\partial [\underline{M}_i]_3}{\partial x_2} \\ 
\frac{\partial [\underline{M}_i]_1}{\partial x_3} + \frac{\partial [\underline{M}_i]_3}{\partial x_1} \\ 
\frac{\partial [\underline{M}_i]_1}{\partial x_2} + \frac{\partial [\underline{M}_i]_2}{\partial x_1} 
\end{bmatrix}
, \quad
\underbrace{\underline{B}_i}_{2\text{D}} = \begin{bmatrix} 
\frac{\partial [\underline{M}_i]_1}{\partial x_1} \\ 
\frac{\partial [\underline{M}_i]_2}{\partial x_2} \\ 
\frac{\partial [\underline{M}_i]_1}{\partial x_2} + \frac{\partial [\underline{M}_i]_2}{\partial x_1} 
\end{bmatrix}
```
"""

# ╔═╡ f53dc83a-af49-4abc-90fc-4b5301987a0a
md"""
With these definitions, we see that we have (in 2d as an example)
```math
\underline{B}_i^\mathrm{T} \underline{\sigma} = 
\frac{\partial [\underline{M}_i]_1}{\partial x_1} \sigma_{11}
+ \frac{\partial [\underline{M}_i]_2}{\partial x_2} \sigma_{22} 
+ \left[\frac{\partial [\underline{M}_i]_1}{\partial x_2} + \frac{\partial [\underline{M}_i]_2}{\partial x_1}\right]\sigma_{12} = \frac{\partial \underline{M}_i}{\partial \underline{x}} :\underline{\underline{\sigma}}
```
Using the Voigt notation, we then get the Finite Element form,
```math
\int_\Omega \underline{M}_i^\mathrm{T} \rho \underline{M}_j\ \mathrm{d}\Omega\ \ddot{a}_j + \int_\Omega \underline{B}_i^\mathrm{T} \underline{\sigma}\ \mathrm{d}\Omega = \int_\Gamma \underline{M}_i^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma - \int_\Omega \underline{M}_i^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega
```
Finally, we insert the constitutive relationship, Hooke's law, ``\underline{\sigma} = \underline{\underline{D}}\ \underline{\varepsilon}``, in combination with using that
```math
\underline{\varepsilon} = \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{bmatrix} = \begin{bmatrix} \frac{\partial u_1}{\partial x_1} \\ \frac{\partial u_2}{\partial x_2}  \\ \frac{\partial u_1}{\partial x_2} + \frac{\partial u_2}{\partial x_1} \end{bmatrix} 
= \sum_{i = 1}^{N_s} \begin{bmatrix} \frac{\partial \left[\underline{M}_i\right]_1}{\partial x_1} \\ \frac{\partial \left[\underline{M}_i\right]_2}{\partial x_2}  \\ \frac{\partial \left[\underline{M}_i\right]_1}{\partial x_2} + \frac{\partial \left[\underline{M}_i\right]_2}{\partial x_1} \end{bmatrix} a_i = \sum_{i = 1}^{N_s} \underline{B}_i a_i = \underline{B}_i a_i
```
we get
```math
\begin{align}
\underbrace{\int_\Omega \underline{M}_i^\mathrm{T} \rho \underline{M}_j\ \mathrm{d}\Omega}_{M_{ij}}\ \ddot{a}_j + \underbrace{\int_\Omega \underline{B}_i^\mathrm{T} \underline{\underline{D}}\ \underline{B}_j\ \mathrm{d}\Omega}_{K_{ij}}\ a_j &= \underbrace{\int_\Gamma \underline{M}_i^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma + \int_\Omega \underline{M}_i^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega}_{f_i} \\ 
M_{ij} \ddot{a}_j + K_{ij} a_j &= f_i \\ 
\underline{\underline{M}}\ \ddot{\underline{a}} + \underline{\underline{K}}\ \underline{a} &= \underline{f}
\end{align}
```
In this course, we will only consider quasi-static problems, i.e. when ``\ddot{\underline{a}} = \underline{0}``, such that our problem becomes
```math
K_{ij} a_j = f_i, \quad \text{or}\quad 
\underline{\underline{K}}\ \underline{a} = \underline{f}
```
!!! note "Vector-valued shape vs mass matrix"
	Please do not mix the mass matrix, ``M_{ij}``, with the vector-valued shape functions ``\underline{M}_i(\underline{x})``. These are different things and not directly related! 

"""

# ╔═╡ cf0b8288-82b5-445a-aa27-76fbdeef35d8
md"""
### Construction of vector-valued shape functions
The vector-valued shape functions, ``\underline{M}_i``, are defined based on our already defined scalar shape functions, ``N_i``, as 
```math
\underline{M}_{2i-1} = \begin{bmatrix} N_i(\underline{x}) \\ 0 \end{bmatrix}, \quad 
\underline{M}_{2i} = \begin{bmatrix} 0 \\ N_i(\underline{x}) \end{bmatrix}
```
And consequently, we get the gradients of these as
```math
\frac{\partial \underline{M}_{2i-1}}{\partial \underline{x}} = \begin{bmatrix} \frac{\partial N_i}{\partial x_1} & 0 \\ 
\frac{\partial N_i}{\partial x_2} & 0 \end{bmatrix}, \quad 
\frac{\partial \underline{M}_{2i}}{\partial \underline{x}} = \begin{bmatrix} 
0 & \frac{\partial N_i}{\partial x_1} \\ 
0 & \frac{\partial N_i}{\partial x_2} \end{bmatrix}
```
And then using the definition of the ``\underline{B}_i`` vector, we get
```math
\underline{B}_{2i-1} = \begin{bmatrix}
\frac{\partial N_i}{\partial x_1} \\ 0 \\ 
\frac{\partial N_i}{\partial x_2} \end{bmatrix}, \quad 
\underline{B}_{2i} = \begin{bmatrix}
0 \\ \frac{\partial N_i}{\partial x_2} \\ 
\frac{\partial N_i}{\partial x_1} \end{bmatrix}
```
You will implement this construction of vector-valued shape functions and their gradients as
* `vectorize_shape_values`
* `voigtize_shape_gradients`

For the latter, the function `tovoigt_strain`, which converts a matrix, ``\partial \underline{M}_i / \partial \underline{x}`` to the voigt vector, ``\underline{B}_i``, can be useful! As usual, we will do all of this on the element numbering level, i.e. we work with ``\underline{B}_i^e`` and ``\underline{M}_i^e``. Furthermore, note that all the parametric mapping etc. is done before doing this construction, i.e. the workflow for a given quadrature point ``\underline{\xi}_q`` is,
1. Calculate the scalar shape values, ``\hat{N}_i^e(\underline{\xi}_q)``
2. Calculate the reference shape gradient values, ``\partial \hat{N}_i^e / \partial \underline{\xi}`` at ``\underline{\xi} = \underline{\xi}_q``. 
3. Calculate the jacobian, ``\underline{\underline{J}} = \partial \underline{x}/\partial \underline{\xi}``
4. Calculate the the total integration weight, ``w_q * \det(\underline{\underline{J}})``
5. Map the gradients, ``\partial N^e_i / \partial \underline{x} = \underline{\underline{J}}^{-T} \partial \hat{N}_i^e / \partial \underline{\xi}``
6. Construct the vectorized shape values, ``\underline{M}^e`` from ``N^e(\underline{\xi})`` using `vectorize_shape_values`
7. Construct the voigtized shape gradients, ``\underline{B}^e``, from ``\partial N^e/\partial \underline{x}`` using `voigtize_shape_gradients`
8. Calculate the contributions to the FE form (or similar for e.g. postprocessing)

"""

# ╔═╡ 72c42622-fa57-4020-aa51-77fd3bd9f79c
md"""
### Global and local numbering
For scalar problems, we chose to use the global numbering scheme that the Degree of Freedom (DoF) number is the same as the node number. For vector-valued problems such as the mechanical equilibrium, where we look for the vector-valued displacement field, ``\underline{u}(\underline{x})``, we have `dim` (e.g. 2 in 2d) DoFs per node, and we cannot simply say DoF number = node number. So for an element with `num_element_nodes` nodes, we have `num_element_dofs = 2 * num_element_nodes` DoFs. As a convention, we get our element dofs, `edofs`, for an element with node number `enods` as
```
edofs_x = 2 * enods - 1;
edofs_y = 2 * enods;
edofs(1:2:num_element_dofs) = edofs_x;
edofs(2:2:num_element_dofs) = edofs_y;
```

!!! note "Alternative convention"
	The chosen convention is just one of many possible. An alternative, given `num_nodes` nodes in the mesh, that we will **not** use, is 
	```
	edofs_x = enods;
	edofs_y = enods + num_nodes;
	edofs(1:2:num_element_dofs) = edofs_x;
	edofs(2:2:num_element_dofs) = edofs_y;
	```
	This will affect how we can translate from knowing the node number to knowing the DoF number of each component. For the remaining of the course, we will **not** use the `edofs_x = enods` and `edofs_y = enods + num_nodes` convention, but the `edofs_x = 2 * enods - 1` and `edofs_y = 2 * enods` convention!

Using the convention `edofs_x = 2 * enods - 1` and `edofs_y = 2 * enods`, we can always get the DoF numbers for the x and y displacements from the node number, `nodenr` as 
```
xdof = 2 * nodenr - 1;
ydof = 2 * nodenr;
```
"""

# ╔═╡ 15b6bf2a-88a4-4782-9087-349a5e47f75b
md"""
### Dirichlet BC
For the chosen numbering convention, we have the properties of shape function nr ``i``, ``\underline{N}_i(\underline{x})``, at node number ``j`` with coordinates ``\underline{x}_j``,
```math
\underline{M}_i(\underline{x}_j) = \left\lbrace \begin{matrix}
[1, 0]^\mathrm{T}, & i = 2j - 1 \\
[0, 1]^\mathrm{T}, & i = 2j\ \ \phantom{-1} \\
[0, 0]^\mathrm{T}, & \text{else}\phantom{--..} % Ugly spacing hack...
\end{matrix} \right.
```

This implies, that when we approximate a function as
```math
\underline{u}(\underline{x}) \approx \underline{u}_h(\underline{x}) = \sum_{i = 1}^{N_\mathrm{dofs}} \underline{M}_i(\underline{x}) a_i
```
we have that 
```math
\underline{u}_h(\underline{x}_j) = \begin{bmatrix} a_{2j-1} \\ a_{2j} \end{bmatrix}
```
So if want to prescribe the x-displacement components at the nodes `dbc_x_nodes` to `ux_c`, we can simply do
```
cdofs = 2 * dbc_x_nodes - 1;
ac = ones(length(dbc_x_nodes), 1) * ux_c;
```
If we also want to prescribe the y displacements at the (possibly) different nodes `dbc_y_nodes` to `uy_c`, we can extend these vectors by
```
cdofs = [cdofs; 2 * dbc_y_nodes];
ac = [ac; ones(length(dbc_y_nodes)) * uy_c];
```
If we have further constraints (e.g. with different constrained values), we can add them by further extending the vector. By knowing the node number, we can even insert constrained values that depend on the coordinates, e.g. let's say that we want to prescribe that ``u_1(\underline{x}) = x_2/10`` on the boundary ``\Gamma_\mathrm{right}`` with nodes `left_nodes`, then we can add this as
```
ux_left = node_coordinates(2, left_nodes);
cdofs = [cdofs; 2 * left_nodes - 1];
ac = [ac; ux_left]
```
Once we have built up all the constrained dofs into the vector `cdofs`, and the corresponding constrained values into `ac`, we can solve for the unknown dof-values, `af`, by doing,
```
fdofs = setdiff((1:ndofs)', cdofs)
af = K(fdofs, fdofs) \ (f(fdofs) - K(fdofs, cdofs) * ac);
```
just as for scalar problems. 
"""

# ╔═╡ 23d40f9d-c661-4551-90a9-e815f9edfca6
md"""
### Neumann BC
Just as for scalar problems, including the Neumann boundary condition requires calculating the contribution,
```math
f_i^\mathrm{NBC} = \int_{\Gamma_\mathrm{NBC}} \underline{M}_i^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma
```
to the global load vector. Similar to the scalar case, we can split the integrals into known and unknown parts on the boundary, i.e.
```math
\int_{\Gamma} \underline{M}_i^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma = \int_{\Gamma_\mathrm{NBC}} \underline{M}_i^\mathrm{T}\ \underline{t}_\mathrm{NBC}\ \mathrm{d}\Gamma + \int_{\Gamma_\mathrm{DBC}} \underline{M}_i^\mathrm{T}\ \underline{t}_\mathrm{DBC}\ \mathrm{d}\Gamma
```
where ``\underline{t}_\mathrm{NBC}`` is known and ``\underline{t}_\mathrm{DBC}`` is unknown. Here we make the split just based on the geometry, but in practice we can also have boundaries where we only know the ``x``-component of the traction but not the ``y``-component (i.e. we know the ``y``-component of the displacement). In practice, that is not a problem, since in that case we just set ``t_2 = 0`` when adding the Neumann BC, and solve using the ``y``-dofs as constrained. In the end, the application of Neumann BC thus requires us to evaluate the integral, 
```math
f_i^\mathrm{NBC} = \int_{\Gamma_\mathrm{NBC}} \underline{M}_i^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma
```
and we will do this by taking the sum over all facets, just as for the scalar heat problem. The only difference introduced in the implementation is that the traction is now a vector instead of the scalar ``q_\mathrm{n}``. For convenience, we also define the `neumann_routine`s such that we can pass both a full traction, ``\underline{t}_\mathrm{f}``, a normal traction, ``t_\mathrm{n}``, such that we apply the total traction
```math
\underline{t} = \underline{t}_\mathrm{f} + t_\mathrm{n}\ \underline{n}
```
by using the normalized normal vector, ``\underline{n} = \underline{n}_\mathrm{w} / ||\underline{n}_\mathrm{n}||``. This simplifies applying e.g. a pressure load, where the traction is normal to the surface. Having calculated this, we can now simply use the same procedure as for scalar problems to approximate the integral numerically, although with the same adaptations as for the assembly of an element (i.e. calculating the vectorized shape values based on the scalar ones).
"""

# ╔═╡ c144eb46-727a-428c-b2a1-f7ce63f2be7a
md"""
### Symmetry BC
If **both** geometry and load is symmetric, we can apply symmetry boundary conditions. That is, the displacement perpendicular to the symmetry line is zero (i.e. zero Dirichlet) and there is zero shear tractions on the symmetry line (i.e. zero Neumann). In our `MATLAB` implementations, we will always have that the symmetry lines are aligned with the global coordinate system, such that we can apply these constraints to each component of the node displacements.

To apply symmetry (or other Dirichlet) boundary conditions in different coordinate systems, we have to use so-called Affine constraints, which creates a coupling between degree's of freedom, but doing this in `MATLAB` is out of scope for this course. However, in the `Abaqus` implementation, you will use a local coordinate system to apply symmetry conditions on a slanted face. 
"""

# ╔═╡ f61e8d2d-1842-4ba8-aff8-01f6b9ad9025
md"""
## L15a: Postprocessing
Solving the boundary value problem implies knowing the dofs, ``a_i``, in the approximation ``\underline{u}(\underline{x}) \approx \underline{M}_i(\underline{x})a_i``. After this, we typically want to perform so-called postprocessing. Typically, we want to calculate the stress in our domain, or reaction forces. 

### Calculating the stress
When calculating the stress, we typically do this in the quadrature points. However, we can also calculate the stress at other points. A typical example that we will use in this course, is that independent of how many quadrature points we use when solving the problem, we sometimes just want to calculate the stress at the center of the element. To calculate the stress, we first need the strain. This we can approximate as 
```math
\underline{\varepsilon}(\underline{x}) \approx \underline{B}_i(\underline{x}) a_i
```
However, as usual we will use the local numbering working on the element level, i.e. 
```math
\underline{\varepsilon}(\underline{\xi}) \approx \underline{B}^e_i(\underline{\xi}) a_i^e
```
such that we can calculate this for the quadrature point ``\underline{\xi}_q`` (note that we don't need the weight ``w_q``, as we do not integrate something). Given the strain, we can then calculate the stress using Hooke's law with the suitable stiffness matrix, ``\underline{\underline{D}}``,
```math
\underline{\sigma} = \underline{\underline D}\ \underline{\varepsilon}
```
"""

# ╔═╡ 33b6a299-29cc-44bf-9e9e-cba48d0b427d
md"""
### Calculating reaction forces
Here, we are interested in calculating a component, ``k\in[1,2]``, of the reaction force, ``F_k``, on a part of the boundary, ``\Gamma_\mathrm{R}``, where we have applied Dirichlet BCs for component ``k``. So specifically, we are interested in calculating
```math
F_k = \int_{\Gamma_\mathrm{R}} t_k\ \mathrm{d}\Gamma
```

Similar to the scalar case, after solving the equation systems we know ``f_i`` for all ``i``, and we have
```math
\begin{align}
f_i &= \int_\Omega \underline{M}_i(\underline{x})^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega + 
\int_\Gamma \underline{M}_i(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma \\
 &= \int_\Omega \underline{M}_i(\underline{x})^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega + 
\int_{\Gamma_\mathrm{R}} \underline{M}_i(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma
 + 
\int_{\Gamma_\text{other}} \underline{M}_i(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma
\end{align}
```
where ``\Gamma_\text{other}`` is the remaining of the boundary, i.e. ``\Gamma_\text{other} = \Gamma \setminus \Gamma_\mathrm{R}``. 

As the scalar shape functions fulfill the following property,
```math
\sum_{i = 1}^{N_\mathrm{s}} N_i(\underline{x}) = 1 \quad \forall \underline{x}\in\Omega
```
the vectorized shape functions fullfill
```math
\sum_{i = 1}^{N_\mathrm{nodes}} \underline{M}_{2i-1}(\underline{x}) = \underline{e}_1 = \begin{bmatrix} 1 \\ 0 \end{bmatrix} \quad \forall \underline{x}\in\Omega, \quad 
\sum_{i = 1}^{N_\mathrm{nodes}} \underline{M}_{2i}(\underline{x}) = \underline{e}_2 = \begin{bmatrix} 0 \\ 1 \end{bmatrix} \quad \forall \underline{x}\in\Omega, \quad 
```
for our numbering convention. 

We also have that 
```math
N_i(\underline{x}) = 0\ \forall\ \underline{x}\in\Gamma_\mathrm{R}\quad \text{if}\quad \underline{x}_i \notin \Gamma_\mathrm{R}
```
This is the same as for the reaction flux earlier. For vector-valued problems this leads to 
```math
\underline{M}_{2i-1}(\underline{x}) = \underline{M}_{2i}(\underline{x}) = \underline{0}\ \forall\ \underline{x}\in\Gamma_\mathrm{R}\quad \text{if}\quad \underline{x}_i \notin \Gamma_\mathrm{R}
```

Together, this implies that
```math
\begin{align}
F_1 &= \int_{\Gamma_\mathrm{R}} t_1\ \mathrm{d}\Gamma = \sum_{i \in \mathbb{R}} \int_{\Gamma_\mathrm{R}} \underline{M}_{2i-1}(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma \\
F_2 &= \int_{\Gamma_\mathrm{R}} t_2\ \mathrm{d}\Gamma = \sum_{i \in \mathbb{R}} \int_{\Gamma_\mathrm{R}} \underline{M}_{2i}(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma 
\end{align}
```
Inserting the known force vector, we get
```math
\begin{align}
F_1 &= \sum_{i \in \mathbb{R}} \left[ f_{2i-1} - 
\underbrace{\int_\Omega \underline{M}_{2i-1}(\underline{x})^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega -
\int_{\Gamma_\text{other}} \underline{M}_{2i-1}(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma}_{f_{2i-1}^\mathrm{known}} \right]
\\
F_2 &= \sum_{i \in \mathbb{R}} \left[ f_{2i} - 
\underbrace{\int_\Omega \underline{M}_{2i}(\underline{x})^\mathrm{T}\ \underline{b}\ \mathrm{d}\Omega -
\int_{\Gamma_\text{other}} \underline{M}_{2i}(\underline{x})^\mathrm{T}\ \underline{t}\ \mathrm{d}\Gamma}_{f_{2i}^\mathrm{known}} \right]
\end{align}
```

For mechanical problems, we will skip the complication of neighboring Dirichlet boundaries, and only consider cases when the traction on is known for nodes in ``\mathbb{R}`` on ``\Gamma_\text{other}``. Consequently, ``f_{2i-1}^\mathrm{known}`` and ``f_{2i}^\mathrm{known}`` are known after assembly. The reaction force in the ``x``-direction, `F1`, from the nodes `rnodes_x` where we have constrained the ``x``-displacements, and the reaction force in the ``y``-direction, `F2`, from the nodes `rnodes_y`, where we have constrained the ``y``-displacements, is then calculated as
```
[K, f_known] = calculate_matrix_and_vector(...) % Add all known contributions, 
										   		% i.e. except unknown tractions
a(cdofs) = ac; % Set constrained displacements
a(fdofs) = K(fdofs, fdofs) \ (f_known(fdofs) - K(fdofs, cdofs) * ac);
f(fdofs) = f_known(fdofs)
f(cdofs) = K(cdofs, :) * a;
rxdofs = 2 * rnodes_x - 1; 
F1 = sum(f(rxdofs) - f_known(rxdofs));
rydofs = 2 * rnodes_y; 
F2 = sum(f(rydofs) - f_known(rydofs));
```
"""

# ╔═╡ 85374012-787e-4ac2-8a20-cf584eb232bc
md"""
## L15b: Robin BC
So far, we have only discussed the standard boundary conditions,
* Dirichlet BC: The primary field (e.g. temperature or displacement) is known
* Neumann BC: The boundary flux or traction is known

However, we can also have a mixed case. For thermal simulations, this corresponds to a heat transfer problem, where the boundary heat flux depends on the difference between the temperature in the structure and the surroundings, i.e. we have
```math
q_\mathrm{n} = \alpha_c [T_\mathrm{s} - T(\underline{x})]
```
where the heat convection coefficient, ``\alpha_c``, is a parameter that depends on the air flow (and thus if the interface is horizontal or vertial), ``T_\mathrm{s}`` is the known temperature of the surroundings, and ``T(\underline{x})`` is the temperature in the structure at the boundary. Starting from the FE form for the heat equation,
```math
\int_\Omega \mathrm{grad}(N_i(\underline{x}))^\mathrm{T} \underline{\underline{D}}\ \mathrm{grad}(N_j(\underline{x})) \ \mathrm{d}\Omega\ a_j = 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - 
\int_\Gamma N_i(\underline{x})\ q_\mathrm{n}\ \mathrm{d}\Gamma
```
we split the boundary terms into the Dirichlet, ``\Gamma_\mathrm{DBC}``, Neumann, ``\Gamma_{NBC}``, and Robin, ``\Gamma_\mathrm{RBC}``, parts. This gives,
```math
\begin{align}
&\int_\Omega \mathrm{grad}(N_i(\underline{x}))^\mathrm{T} \underline{\underline{D}}\ \mathrm{grad}(N_j(\underline{x})) \ \mathrm{d}\Omega\ a_j 
= 
\int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega \\ 
&- \int_{\Gamma_\mathrm{DBC}} N_i(\underline{x})\ q_\mathrm{n}^\mathrm{DBC}\ \mathrm{d}\Gamma
- \int_{\Gamma_\mathrm{NBC}} N_i(\underline{x})\ q_\mathrm{n}^\mathrm{NBC}\ \mathrm{d}\Gamma
- \int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c[T_\mathrm{s} - T(\underline{x})]\ \mathrm{d}\Gamma
\end{align}
```
where ``q_\mathrm{n}^\mathrm{NBC}`` is known while ``q_\mathrm{n}^\mathrm{DBC}`` is unknown. If we insert the FE approximation in the Robin BC, we get
```math
\int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c[T_\mathrm{s} - T(\underline{x})]\ \mathrm{d}\Gamma = \int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c T_\mathrm{s}\ \mathrm{d}\Gamma - \underbrace{\int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c N_j(\underline{x})\ \mathrm{d}\Gamma}_{K_{ij}^c}\ a_j
```
and we notice that we get another contribution to the stiffness matrix, ``K_{ij}^c``. If we put all of this together in the complete FE form, we now have
```math
\begin{align}
& \underbrace{\int_\Omega \mathrm{grad}(N_i(\underline{x}))^\mathrm{T} \underline{\underline{D}}\ \mathrm{grad}(N_j(\underline{x})) \ \mathrm{d}\Omega
+ \int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c N_j(\underline{x})\ \mathrm{d}\Gamma}_{K_{ij} + K_{ij}^c}\ a_j \\ 
&= \int_\Omega N_i(\underline{x})\ h\ \mathrm{d}\Omega - \int_{\Gamma_\mathrm{DBC}} N_i(\underline{x})\ q_\mathrm{n}^\mathrm{DBC}\ \mathrm{d}\Gamma
- \int_{\Gamma_\mathrm{NBC}} N_i(\underline{x})\ q_\mathrm{n}^\mathrm{NBC}\ \mathrm{d}\Gamma
- \int_{\Gamma_\mathrm{RBC}} N_i(\underline{x})\ \alpha_c\ T_\mathrm{s}\ \mathrm{d}\Gamma
\end{align}
```

To implement Robin-type boundary conditions, we have to combine our knowledge for how to assemble the stiffness matrix (element routines) with how to integrate on the boundary (neumann routines). But the same applies: we want to calculate the local matrix,
```math
K_{ij}^{c,e} = \int_{\Gamma^f} N^e_i(\underline{x})\ \alpha_c N^e_j(\underline{x})\ \mathrm{d}\Gamma
```
and the local vector,
```math
f_i^e = \int_{\Gamma^f} N^e_i(\underline{x})\ \alpha_c\ T_\mathrm{s}\ \mathrm{d}\Gamma
```
for the facet ``f`` belonging to element ``e``.

"""

# ╔═╡ d8a8fc12-ea0c-4785-a34e-3797aeb6be0c
md"""
We can also introduce Robin BCs for the mechanical problem, this corresponds to an elastic bed, which can be useful for analyzing e.g. a beam lying on the ground. Here, we typically have the traction given as
```math
\underline{t} = -k\ \underline{u}
```
(although it could also be ``\underline{t} = k[\underline{u}_s - \underline{u}]``).
"""

# ╔═╡ 7e7cf116-f844-4345-bb00-e9829e6262be
begin
	"""
		eq(equation_text, label)
	
	Write an numbered text block (typically an equation) which can be referred to with `[^label]` in other parts of the document.
	"""
	function eq(equation_text, label)
		labelstr = "[^$label]:"
		Columns(equation_text, Markdown.parse(labelstr); widths=[100, 0], gap = 0)
	end
end;

# ╔═╡ d8328904-bc46-4c04-a0ff-6bf34bdcdc39
eq(md"""
```math
	\begin{align}
	    x(\xi) = \sum_{\alpha=1}^{N_\mathrm{nodes}} \hat{N}_\alpha(\xi) x_\alpha
	\end{align}
```""",
   "parametricelement1d")

# ╔═╡ fe0acbb2-f60d-43d9-8694-f4d6ef2a0025
eq(md"""
```math
   \mathrm{div}(\underline{q}) = \frac{\partial q_1}{\partial x_1} + \frac{\partial q_2}{\partial x_2} = \frac{\partial q_i}{\partial x_i}
```""", "divergence")

# ╔═╡ 4b7acd47-9515-4325-bfb3-2899f847c26b
eq(md"""
```math
\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_\Omega \mathrm{div}(\underline{q})\ \mathrm{d}\Omega
```""", "divergencetheorem")

# ╔═╡ c1e878bd-0aac-4f28-a528-75c9ce9e3e4a
eq(md"""
```math
\mathrm{div}(\underline{q}) = h
```""", "heatequationstrong")

# ╔═╡ 64a938b2-6001-4424-84e9-c599020fde0b
eq(md"""
   ```math
   \int_\Omega \mathrm{div}(\delta T\ \underline{q})\ \mathrm{d}\Omega = 
   \int_\Gamma \delta T\ \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_\Omega \mathrm{grad}(\delta T)\cdot \underline{q}\ \mathrm{d}\Omega + \int_\Omega \delta T\ \mathrm{div}(\underline{q})\ \mathrm{d}\Omega
   ```
   """, "greengaussthm")

# ╔═╡ 70c347f5-f365-434b-a557-15d84b45997c
eq(md"""
```math
\int_{\Gamma^f} h(\underline{x})\ \mathrm{d}\Gamma \approx \sum_{q=1}^{N_\mathrm{qp}} h(\underline{x}) \left\vert\underline{n}_\mathrm{w}^f\right\vert\ w_q^\mathrm{line}
```
""", "facetintegration")

# ╔═╡ 6bfc5b5e-cb6f-4e33-8b86-f99ffa713eb5
begin
"""
	topoints(x)::Vector{Point{1}}
	topoints(x, y)::Vector{Point{2}}
	topoints(x, y, z)::Vector{Point{3}}

Convert a set of coordinate vectors into a vector of `GeometryBasics.Point`s
"""
function topoints(args::Vararg{Any,N}) where N
    T = promote_type(eltype.(args)...)
    points = Vector{Plt.Point{N, T}}(undef, length(first(args)))
    map!(Plt.Point, points, args...)
    return points
end
"""
	topoints([T = Float64], ::Val{dim})

Quick function to setup an empty set of points for creating the observable
"""
topoints(::Val{dim}) where {dim} = topoints(Float64, Val(dim))
function topoints(::Type{T}, ::Val{dim}) where {T<:Number, dim}
	return topoints(ntuple(_ -> zeros(T, 1), dim)...)
end
function meshgrid(xv, yv)
    X = reshape([x for x in xv, _ in 1:length(yv)], :)
    Y = reshape([y for _ in 1:length(xv), y in yv], :)
    return collect(X), collect(Y)
end
end;

# ╔═╡ df75329d-1ab5-443c-8eea-38060b4461d9
begin
	a1_base2_slider = @bind a1_base2 Slider(-1.2:0.1:1.2; default=0, show_value = true)
	a2_base2_slider = @bind a2_base2 Slider(-1.2:0.1:1.2; default=0, show_value = true)
	function make_base2_plot()
		fig = Plt.Figure(size=(500,250))
		ax = Plt.Axis(fig[1,1]; xlabel = L"\xi", ylabel = L"f(\xi)")
		points = Plt.Observable(topoints(Val(2)))
		Plt.lines!(ax, points)
		Plt.scatter!(ax, points)
		Plt.xlims!(ax, 1.05 .* (-1, 1))
		Plt.ylims!(ax, 1.05 .* (-1.2, 1.2))
		return fig, points
	end
	fig_base2, data_base2 = make_base2_plot()
	md"""
	By calculating the sum, 
	```math
	f(\xi;\underline{a}) = \hat{N}_1(\xi) a_1 + \hat{N}_2(\xi) a_2 = \sum_{i=1}^2 \hat{N}_i(\xi) a_i = \hat{N}_i(\xi) a_i
	```
	we can approximate different functions by adjusting the DoFs, ``a_1`` and ``a_2``.

	``a_1``: $(a1_base2_slider)
	
	``a_2``: $(a2_base2_slider)
	"""
end

# ╔═╡ 9f9e4e47-42ef-4646-bcc2-3d7e4772114c
begin
	data_base2[] = topoints([-1.0, 1.0], [a1_base2, a2_base2])
	fig_base2
end

# ╔═╡ f2f0cafd-e220-4555-8029-47adcbe6d425
begin
	a1_base3_slider = @bind a1_base3 Slider(-1.2:0.1:1.2; default=0, show_value = true)
	a2_base3_slider = @bind a2_base3 Slider(-1.2:0.1:1.2; default=0, show_value = true)
	a3_base3_slider = @bind a3_base3 Slider(-1.2:0.1:1.2; default=0, show_value = true)
	function make_base3_plot()
		fig = Plt.Figure(size=(500,250))
		ax = Plt.Axis(fig[1,1]; xlabel = L"\xi", ylabel = L"f(\xi)")
		points = Plt.Observable(topoints(Val(2)))
		nodes = Plt.Observable(topoints(Val(2)))
		Plt.lines!(ax, points)
		Plt.scatter!(ax, nodes)
		Plt.xlims!(ax, 1.05 .* (-1, 1))
		Plt.ylims!(ax, 1.05 .* (-1.5, 1.5))
		return fig, (; points, nodes)
	end
	fig_base3, data_base3 = make_base3_plot()
	md"""
	By calculating the sum, 
	```math
	f(\xi;\underline{a}) = \hat{N}_1(\xi) a_1 + \hat{N}_2(\xi) a_2 + \hat{N}_3(\xi) a_3 = \sum_{i=1}^3 \hat{N}_i(\xi) a_i = \hat{N}_i(\xi) a_i
	```
	we can approximate different functions by adjusting the DoFs, ``a_1``, ``a_2``, and ``a_3``.

	``a_1``: $(a1_base3_slider)
	
	``a_2``: $(a2_base3_slider)

	``a_3``: $(a3_base3_slider)
	"""
end

# ╔═╡ 4907273f-8cc4-4e41-b5be-ad4ccfc78822
begin
	data_base3.points[] = (ξ = collect(range(-1, 1, 50));
		topoints(ξ, map(ξ) do x
			a = (a1_base3, a2_base3, a3_base3)
			ip = Lagrange{RefLine,2}()
			sum(i -> a[i] * Ferrite.reference_shape_value(ip, Vec{1}((x,)), i), 1:3; init = 0.0)
		end))
	data_base3.nodes[] = topoints([-1.0, 1.0, 0.0], [a1_base3, a2_base3, a3_base3])
	fig_base3
end

# ╔═╡ f436cc3b-49cb-4f7b-915e-e2c0fee4c2cb
FootnotesNumbered()

# ╔═╡ 87f373ca-78eb-4e1f-b539-33d3fd05ab21
FootnotesInlineStyleBaseline()

# ╔═╡ 6444d515-be51-40d4-b088-fdca80a4ec72
begin
	tongonface(c::Union{Triangle, Quadrilateral}) = GB.NgonFace(Ferrite.get_node_ids(c))
	topoint(n::Node{sdim}) where {sdim} = GB.Point{sdim}(get_node_coordinate(n))

	function tomesh(grid::Ferrite.AbstractGrid)
	    faces = map(tongonface, getcells(grid))
	    vertices = map(topoint, getnodes(grid))
	    return GB.Mesh(vertices, faces)
	end
end;

# ╔═╡ cd3d341f-8517-44f9-add8-d844f2e5735a
begin
	function calculate_vonmises(m, _, ∇u, old)
		σ, _, _ = material_response(m, symmetric(∇u), old)
		return sqrt((3/2) * σ ⊡ dev(σ))
	end
		
	function showcase_mech(t)
		grid = generate_grid(Quadrilateral, (40, 10), zero(Vec{2}), Vec((4.0, 1.0))) # TODO: Create inp file with more interesting geometry
		ip = Lagrange{RefQuadrilateral, 2}()^2
		dh = close!(add!(DofHandler(grid), :u, ip))
		m = ReducedStressState(PlaneStrain(), LinearElastic(;E = 80.e3, ν = 0.2))
		qr = QuadratureRule{RefQuadrilateral}(2)
		cv = CellValues(qr, ip)
		db = setup_domainbuffer(DomainSpec(dh, m, cv))
		K = allocate_matrix(dh)
		f = zeros(ndofs(dh))
		ch = ConstraintHandler(dh)
		add!(ch, Dirichlet(:u, getfacetset(grid, "left"), Returns(zero(Vec{2}))))
		add!(ch, Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> t, 2))
		close!(ch)
		work!(start_assemble(K, f), db);
		proj = L2Projector(grid)
		add!(proj, 1:getncells(grid), Lagrange{RefQuadrilateral, 1}(); qr_rhs = qr)
		close!(proj)
		qpeval = QuadPointEvaluator{Float64}(db, calculate_vonmises)
		mesh = tomesh(grid)

		fill!(f, 0)
		update!(ch, t)
		apply!(K, f, ch)
		a = K \ f;
		work!(qpeval, db; a)
		σvm = project(proj, qpeval.data)
		σvm_nodes = evaluate_at_grid_nodes(proj, σvm)
		u_nodes = evaluate_at_grid_nodes(dh, a, :u)
		map!(mesh.position, mesh.position, u_nodes) do X, u
			X + GB.Vec2(u)
		end
		return mesh, σvm_nodes
	end
	function showcase_mech_setup_plot()
		fig = Plt.Figure(size = (800, 400))
		ax = Plt.Axis(fig[1,1]; aspect = Plt.DataAspect(), xlabel = "x [mm]", ylabel = "y [mm]")
		Plt.xlims!(ax, -0.25, 4.25)
		Plt.ylims!(ax, -0.6, 1.6)
		mesh, σvm_nodes = showcase_mech(0.0)
		
		mesh_plt = Plt.mesh!(ax, mesh, color = σvm_nodes, colorrange = (0, 5000))
		wire_plt = Plt.wireframe!(ax, mesh; color = :black, linewidth = 1)
		Plt.Colorbar(fig[1,2], mesh_plt; colorrange = (0, 5000), label = "von Mises stress [MPa]")
		return fig, (; mesh_plt, wire_plt)
	end
	function showcase_mech_update_plot!(pltdata; t)
		mesh, σvm_nodes = showcase_mech(t)
		Plt.update!(pltdata.mesh_plt; arg1 = mesh, color = σvm_nodes)
		Plt.update!(pltdata.wire_plt; arg1 = mesh)
	end
	showcase_mech_fig, showcase_mech_data = showcase_mech_setup_plot()
	mech_showcase_load_slider = @bind mech_showcase_load Slider(-0.5:0.1:0.5; default = 0, show_value = true)
	
	md"""
	## L12a: Mechanical equilibrium
	So far, we have considered only the heat equation, now we will apply the same theory to mechanical equilibrium. To do so, we start by defining the stress according to *Cauchy's stress theorem* before determining the force and torque equilibrium conditions. In the 2nd half, L12b, we introduce linear elasticity for both 3d and 2d cases.
	
	But first, let's look at the result of a mechanical simulation with linear elasticity, and you can explore how the boundary conditions affect the effective (von Mises) stress response:
	
	Prescribed ``u_2`` on right boundary: $(mech_showcase_load_slider)
	"""
end

# ╔═╡ ef6fe929-9305-4a9a-84c4-cfbd601a91c4
begin
	showcase_mech_update_plot!(showcase_mech_data; t = mech_showcase_load)
	showcase_mech_fig
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Ferrite = "c061ca5d-56c9-439f-9c0e-210fe06d3992"
FerriteAssembly = "fd21fc07-c509-4fe1-9468-19963fd5935d"
Format = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
MaterialModelsBase = "af893363-701d-44dc-8b1e-d9a2c129bfc9"
MechanicalMaterialModels = "b3282f9b-607f-4337-ab95-e5488ab5652c"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
CairoMakie = "~0.15.7"
Ferrite = "~1.2.0"
FerriteAssembly = "~0.3.7"
Format = "~1.3.7"
GeometryBasics = "~0.5.10"
LaTeXStrings = "~1.4.0"
MaterialModelsBase = "~0.3.1"
MechanicalMaterialModels = "~0.3.0"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "9ee2a05719d0613131abccfdc06562be9a018858"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "1778fd03576b0b6f88d0eafe89c54a3fb8df96a3"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "21f3ae106d1dcc20a66e96366012f7289ebba498"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.5"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "783b21581a051ac91a3921ee37e26a23ed7f57a6"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.5"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.Ferrite]]
deps = ["EnumX", "ForwardDiff", "LinearAlgebra", "NearestNeighbors", "OrderedCollections", "Preferences", "Reexport", "SparseArrays", "StaticArrays", "Tensors", "WriteVTK"]
git-tree-sha1 = "26e13c8e1ed9df75129be2e22409e93ba5565ce1"
uuid = "c061ca5d-56c9-439f-9c0e-210fe06d3992"
version = "1.2.0"

    [deps.Ferrite.extensions]
    FerriteBlockArrays = "BlockArrays"
    FerriteMetis = "Metis"
    FerriteSparseMatrixCSR = "SparseMatricesCSR"

    [deps.Ferrite.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b"
    SparseMatricesCSR = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"

[[deps.FerriteAssembly]]
deps = ["ConstructionBase", "Ferrite", "ForwardDiff", "MaterialModelsBase"]
git-tree-sha1 = "755be5aa347d806f665ee5a076258baba0b56ea8"
uuid = "fd21fc07-c509-4fe1-9468-19963fd5935d"
version = "0.3.7"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "d60eb76f37d7e5a40cc2e7c36974d864b82dc802"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.1"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bfcd42851cf2f1b303f51525a54dc5e98d408a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.15.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cd33c7538e68650bd0ddbb3f5bd50a4a0fa95b50"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.0"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "02b61501dbe6da3b927cc25dacd7ce32390ee970"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.2"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "aa971a09f0f1fe92fe772713a564aa48abe510df"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.3"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "7e6151c8432b91e76d9f9bc3adc6bbaecd00ec0a"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.7"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaterialModelsBase]]
deps = ["ForwardDiff", "StaticArrays", "Tensors"]
git-tree-sha1 = "0f53bc668749dce2b5e79a0b38a46bd2ebf2f4e9"
uuid = "af893363-701d-44dc-8b1e-d9a2c129bfc9"
version = "0.3.1"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.MechanicalMaterialModels]]
deps = ["ForwardDiff", "LinearAlgebra", "MaterialModelsBase", "Newton", "StaticArrays", "Tensors"]
git-tree-sha1 = "b224b75af462aa99995dbe0d118d85063cbbf347"
uuid = "b3282f9b-607f-4337-ab95-e5488ab5652c"
version = "0.3.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "e45bb6034fdef63d0c49b82ba9b889215bf8b344"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.24"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Newton]]
deps = ["DiffResults", "ForwardDiff", "LinearAlgebra", "Preferences", "Printf", "StaticArrays", "Tensors"]
git-tree-sha1 = "ff2f9e38116d816f44706c0c47686eaa00f8877a"
uuid = "83aa5b51-0588-403c-85e4-434ec185aae7"
version = "0.2.2"

    [deps.Newton.extensions]
    RecursiveFactorizationExt = "RecursiveFactorization"

    [deps.Newton.weakdeps]
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d922b4d80d1e12c658da7785e754f4796cc1d60d"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.36"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "064b532283c97daae49e544bb9cb413c26511f8c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.8"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Tensors]]
deps = ["ForwardDiff", "LinearAlgebra", "PrecompileTools", "SIMD", "StaticArrays", "Statistics"]
git-tree-sha1 = "399ec4c46786c47380121f8a6497c65b89f0573f"
uuid = "48a634ad-e948-5137-8d70-aa71f2a747f4"
version = "1.16.2"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "83360bda12f61c250835830cc40b64f487cc2230"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.1"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "a329e0b6310244173690d6a4dfc6d1141f9b9370"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.21.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "5c959b708667b34cb758e8d7c6f8e69b94c32deb"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.15.1+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─6b0414f5-13f1-4e68-944d-d29b18affaa3
# ╟─13c0d238-e1cc-473c-a41f-dfbb05622b0f
# ╟─b92cd124-af25-11f0-b55d-45d572d7ae78
# ╟─49a9b55e-16ad-4d80-ac5d-cbf18cf717d7
# ╟─9bda1921-1001-4bcb-9d53-101fe7010b72
# ╟─df75329d-1ab5-443c-8eea-38060b4461d9
# ╟─9f9e4e47-42ef-4646-bcc2-3d7e4772114c
# ╟─42924ddb-99ac-4d9d-8092-8b43b8f6b9f4
# ╟─f2f0cafd-e220-4555-8029-47adcbe6d425
# ╟─4907273f-8cc4-4e41-b5be-ad4ccfc78822
# ╟─6a1731e3-c2bd-41af-8f7b-cc515d2e7912
# ╟─609a415e-0734-4061-997a-4cfb8b71a12f
# ╟─216d9ced-4729-40a0-9a45-702b85dd9c3d
# ╟─88e934f1-2a7f-478c-8876-f7edbfa1088b
# ╟─1051d025-50fe-4311-94a7-606c46a9e84d
# ╟─889968e5-2713-4f0d-9b7c-b83a9180f2ae
# ╟─e134610b-dd57-4f67-bc00-b5207288d381
# ╟─4e1b0052-2f20-453e-a259-b3e9f1b62752
# ╟─d8328904-bc46-4c04-a0ff-6bf34bdcdc39
# ╟─40615160-0e12-4a1c-877a-507d1c7b834b
# ╟─abf84c67-87f2-4b94-b5a3-578b854e0711
# ╟─d55d3180-eb0b-4d8b-afb2-4daa3dee7b65
# ╟─c3d23889-ff92-476d-92d0-726812079777
# ╟─6fd9e877-66dc-4a59-8fe2-6d6c9c5af473
# ╟─da8d686b-ccdf-49b0-9fa0-20b2c5b339c8
# ╟─8b268765-0e20-49fc-b5f8-e12f6f6136b5
# ╟─9e31346f-75ce-4d57-9a3a-644e14d5b7bf
# ╟─a1df8715-af64-4701-9b1e-b5ae42e876e4
# ╟─082d93db-9088-4a41-bfc4-d1578622a76c
# ╟─dc96f2fd-8afc-40b6-afb0-0a613023f7b7
# ╟─6e4845f2-a984-43f0-9ef6-39f868ed9e20
# ╟─8fffc673-54a3-4c5a-bcae-e4af09cab9bf
# ╟─6e02d8e9-faef-4f90-af0a-fc3f091b55fc
# ╟─e60e838b-2f53-45d9-b644-615f7171aa0b
# ╟─7ee86550-b9cd-4f6e-b039-c59d9c50f563
# ╟─177135bc-71b8-4298-9c21-8c8a11c0c753
# ╟─a004f4bc-1a6a-44f5-810b-58d446f813f4
# ╟─3a99bad7-8447-435d-ae76-392361ca656d
# ╟─9b75b720-0924-4933-9e60-3a0662dfeea7
# ╟─de47e61e-f6db-42b2-96a7-8cd8f31d8825
# ╟─aec4b7af-d986-4d33-9f96-dd6c749df0d3
# ╟─0363492e-8ecf-4f6c-ab50-2b0612572d39
# ╟─080f125a-8003-41e7-b50d-05a4de587301
# ╟─58c1a7d8-03fb-409d-8589-46d2fa597f82
# ╟─8786dc4e-f7f4-4ccc-85f5-5302988c696c
# ╟─fe0acbb2-f60d-43d9-8694-f4d6ef2a0025
# ╟─5608fe81-9b47-4214-a8ad-0960ca147b27
# ╟─340f8448-6db6-404f-a7ea-192aa339c286
# ╟─d9cdd761-bcb3-461a-b473-ab1e6493a80b
# ╟─33618c32-a221-487c-b1ac-0b4beb347c7d
# ╟─4b7acd47-9515-4325-bfb3-2899f847c26b
# ╟─ac197737-6f54-4f22-a049-0e268c7ec4db
# ╟─b1999d4a-3111-4203-ab1c-173905f06d80
# ╟─c1e878bd-0aac-4f28-a528-75c9ce9e3e4a
# ╟─33c4264d-78bb-4c95-9d27-aa340df03808
# ╟─fccdcef3-3768-47d4-ba0e-ecca7308d356
# ╟─64a938b2-6001-4424-84e9-c599020fde0b
# ╟─02ebd323-03e2-41cb-903f-cb0bf690cfd5
# ╟─27087e5b-08fb-4e0e-962c-30999860e238
# ╟─d9385605-737f-401b-aba9-5d8d7a18a18c
# ╟─a658416a-108e-4f2f-9cb1-324e5f2960f7
# ╟─32c1a2ed-ea2a-469b-936f-1a2b4b3eec90
# ╟─ef71b38b-ec86-41eb-86f2-ea8087038502
# ╟─986c99b7-81b8-4ceb-9eee-14e2d0d26b34
# ╟─d40d0648-ea58-47d1-b2b7-77c30005f50a
# ╟─e2f75d06-a29f-4709-a3ca-3dc98cd4f570
# ╟─ce1d0468-b97a-4ae1-b93c-28cc10acb7d4
# ╟─76a2f37c-2672-431a-9a33-128681607393
# ╟─903164a8-30d5-466f-b399-eb380a895117
# ╟─dfe56d4b-13c2-4ad4-ab5e-578f470e4589
# ╟─e85161d4-4b28-47d4-a05e-a3c3e34415b5
# ╟─fdc596a9-bd9e-45bd-a082-0418bb5b5084
# ╟─9aec05ec-d7eb-4ef1-a56f-96a4c2371897
# ╟─67c6438f-2547-4580-82fe-a3263203e939
# ╟─a4d60a06-d262-4cd1-aad2-f19dd221b4df
# ╟─8cc38e33-13f2-47a9-8b33-66ce560fbb3b
# ╟─cad4715c-72ef-452c-bbf4-7bbd6149ef13
# ╟─70c347f5-f365-434b-a557-15d84b45997c
# ╟─510466ef-7197-4e08-b915-1a23eb5fd095
# ╟─87c10954-05f0-4866-ae45-df3fcc103c01
# ╟─e840274e-902e-4b28-a8bd-96af78380a20
# ╟─c21b70de-f36d-4c91-9a29-8204bb0f8cc2
# ╟─a471df59-3c75-4cfa-8f19-8612b20c59b0
# ╟─11bd2c69-1744-4789-8d0b-62c23ae7b896
# ╟─55570b97-a717-4a7c-89c4-a2ae941a8e32
# ╟─af3043ce-2863-44d1-8e7e-db5465651c52
# ╟─ef099aa1-dd65-468e-af83-632ef4e0b535
# ╟─509c3307-b814-4da9-aecc-7adb47f8ec95
# ╟─cd3d341f-8517-44f9-add8-d844f2e5735a
# ╟─ef6fe929-9305-4a9a-84c4-cfbd601a91c4
# ╠═1f5672bc-28af-4ab7-b159-197ebf1e12a3
# ╟─4167e37d-b9e2-484f-935c-729e0507630b
# ╟─bebe6b31-f5f3-4719-b970-7f0583bc3674
# ╟─0d31dfd6-df4b-4d08-a099-c69c2f2cb071
# ╟─9459013d-e12c-4a82-9636-5ce912cccaf3
# ╟─48ebeebf-8507-4941-9361-4a5b4fe601db
# ╟─a9e3d9df-5ab5-4ced-a4e6-c61ce2e79046
# ╟─75a525f7-e727-4c55-9fcd-a4335d97c432
# ╟─dc02ad0d-058b-4369-be78-44039aa6e6b1
# ╟─3345a144-6871-4b5b-b60b-db702770b30a
# ╟─b5815863-6379-4711-a9fd-4cc69e8388e1
# ╟─2ec78a16-f64e-4795-8305-d1a36899eb0a
# ╟─6ccce575-56b5-4fc6-969e-b71348c44a81
# ╟─20f8718e-60ea-415b-aab3-17a90c7f31c6
# ╟─d7e7d8c5-0673-4c78-a484-4ee8116f1a76
# ╟─7cec58dd-dd80-4ce5-b885-f1e2527e481d
# ╟─f13237a3-1bf4-4ec9-a797-20458df04e52
# ╟─e2f98607-c9a6-4a89-ba81-87b5a83cb789
# ╟─5f0cf5a1-bd8c-4637-a28d-e4fd134254ab
# ╟─976fbb59-1adf-4e95-84b8-75c7bf302917
# ╟─f53dc83a-af49-4abc-90fc-4b5301987a0a
# ╟─cf0b8288-82b5-445a-aa27-76fbdeef35d8
# ╟─72c42622-fa57-4020-aa51-77fd3bd9f79c
# ╟─15b6bf2a-88a4-4782-9087-349a5e47f75b
# ╟─23d40f9d-c661-4551-90a9-e815f9edfca6
# ╟─c144eb46-727a-428c-b2a1-f7ce63f2be7a
# ╟─f61e8d2d-1842-4ba8-aff8-01f6b9ad9025
# ╟─33b6a299-29cc-44bf-9e9e-cba48d0b427d
# ╟─85374012-787e-4ac2-8a20-cf584eb232bc
# ╟─d8a8fc12-ea0c-4785-a34e-3797aeb6be0c
# ╟─7e7cf116-f844-4345-bb00-e9829e6262be
# ╟─6bfc5b5e-cb6f-4e33-8b86-f99ffa713eb5
# ╟─f436cc3b-49cb-4f7b-915e-e2c0fee4c2cb
# ╟─87f373ca-78eb-4e1f-b539-33d3fd05ab21
# ╟─6444d515-be51-40d4-b088-fdca80a4ec72
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
