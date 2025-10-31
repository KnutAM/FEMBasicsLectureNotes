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
	using PlutoUI
	using PlutoTeachingTools
	using Ferrite
	using LaTeXStrings
	using Printf
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
This material is provided under the MIT License: [license file](https://raw.githubusercontent.com/KnutAM/FEMBasicsLectureNotes/refs/heads/main/LICENSE)
"""

# ╔═╡ 49a9b55e-16ad-4d80-ac5d-cbf18cf717d7
md"""
## L1a: Approximating functions with shape functions

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
    K_{ij} = \left[\int_{-1}^{1} \hat{N}_i(\xi) \hat{N}_j(\xi)\ \mathrm{d}\xi\right], \quad 
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
	## L1b: Numerical integration (1D)
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

Equipped with the ability to integrate any function numerically, we can now calculate the integrals in L1a, to establish ``K_{ij}`` and ``f_i``. Specifically, if we use the `linear_line_reference_shape_values` or `quadratic_line_reference_shape_values` functions in `BasicFEM`, we obtain for a given coordinate ``\xi``, the vector ``\underline{N} = [\hat{N}_1(\xi), \cdots, \hat{N}_{N_\mathrm{s}}(\xi)]``, such that our "function" to integrate to obtain the matrix ``K_{ij}`` becomes `N' * N` in `MATLAB`.

Solving this problem is the first homework assignment
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
	## L2: Introduction to MATLAB
	See notes on Canvas
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
		fig = Plt.Figure()
		ax_ref = Plt.Axis(fig[1,1]; xlabel = L"\xi", ylabel = L"\hat{N}(\xi)")
		Plt.lines!(ax_ref, [-1, 1], [1, 0]; label = L"\hat{N}_1(\xi) = [1-\xi]/2,\ \partial\hat{N}_1/\partial\xi = -1/2")
		Plt.lines!(ax_ref, [-1, 1], [0, 1]; label = L"\hat{N}_2(\xi) = [1+\xi]/2,\ \partial\hat{N}_2/\partial\xi = 1/2")
		Plt.ylims!(ax_ref, -0.1, 1.8)
		Plt.axislegend(ax_ref; position = :ct)
		ax_glob = Plt.Axis(fig[2,1]; xlabel = L"x", ylabel = L"N(x)")
		Plt.xlims!(ax_glob, 0.0, 2.0)
		Plt.ylims!(ax_glob, -0.1, 1.8)
		N1_line = Plt.lines!(ax_glob, [0.0, 1.0], [1.0, 0.0]; label = L"N_1(x),\ \partial N_1/\partial x = -1/L_e")
		N2_line = Plt.lines!(ax_glob, [0.0, 1.0], [0.0, 1.0]; label = L"N_2(x),\ \partial N_2/\partial x = 1/L_e")
		Plt.axislegend(ax_glob; position = :ct)
		return fig, (;N1_line, N2_line)
	end
	parametric_line_gradient_fig, parametric_line_gradient_data = setup_parametric_line_gradient_fig();
	parametric_line_length_slider = @bind parametric_line_length Slider(0.1:0.1:2.0; default = 1, show_value = true)
md"""
### Mapping of gradients
When solving actual Partial Differential Equations (PDEs), we will also need to get the gradient of a function, ``f(x)``, that we approximate using shape functions, ``g(x) \approx f(x; \underline{a}) = \sum_{i=1}^{N_\mathrm{s}} \hat{N}_i(\xi(x)) a_i``. Specifically, we want to calculate
```math
    \frac{\partial g}{\partial x} \approx \frac{\partial f}{\partial x} = \sum_{i = 1}^{N_\mathrm{s}} \frac{\partial \hat{N}_i(\xi(x))}{\partial x} a_i
```
Since our shape functions, ``\hat{N}(\xi)``, are described as a function of the reference coordinate, ``\xi``, we now need the inverse mapping: ``\xi = \xi(x)`` if we want to calculate the physical coordinate dependent shape function, ``N_i(x) = \hat{N}_i(\xi(x))`` (Notice that we remove the hat, ``\hat{\bullet}``, when denoting the shape function described by the physical coordinates). However, we don't have an explicit function for this inverse mapping, and instead we use
```math
    \underline{\underline{I}} = \frac{\partial x}{\partial x} = \frac{\partial x}{\partial \xi}\frac{\partial \xi}{\partial x} \Rightarrow \frac{\partial \xi}{\partial x} = \left[ \frac{\partial x}{\partial \xi} \right]^{-1} = J^{-1} 
```
to obtain the derivative, ``\partial \xi/\partial x``. The gradient of a shape function can be calculated as,
```math
    \frac{\partial N_i}{\partial x} = \frac{\partial \hat{N}_i}{\partial \xi} \frac{\partial \xi}{\partial x} = \frac{\partial \hat{N}_i}{\partial \xi} J^{-1}
```
This means that we can now calculate the reference shape gradients, ``\partial N_i/ \partial\xi``, and use those to obtain the spatial gradients when required.

For a **linear line element**, the jacobian, ``J``, is given by
```math
J = \frac{\partial \hat{N}_1}{\partial\xi}x_1 + \frac{\partial \hat{N}_1}{\partial\xi}x_2 = \frac{-1}{2}x_1 + \frac{1}{2}x_2 = \frac{1}{2}[x_2 - x_1] = \frac{L_e}{2}
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


The second is Fourier's law for thermal transport,
```math
q = -k \frac{\partial T}{\partial x}
```
stating that the heat flux (energy / (time and surface area)) is proportional (with material constant ``k``) to the negative temperature gradient. This implies that energy goes from hot areas to cold areas, as we can observe e.g. when cooking.

### Strong form
To derive the strong form, we need to put the 1st law of thermodynamics into a formula. To do this, we must first define the following quantities

* Internal energy, ``e``: The stored internal energy depends on the temperature, and in the linear case we have that ``e = \rho c_\mathrm{p} T``, where ``\rho`` is the density (mass/volume) and ``c_\mathrm{p}`` the heat capacity (energy / (mass and temperature)). 

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
\begin{bmatrix} a_f \\ a_c \end{bmatrix} = \begin{bmatrix} f_f \\ f_c \end{bmatrix}
```
where we in `MATLAB` can obtain these submatrices and vectors as 
```
K_ff = K(fdofs, fdofs); K_fc = K(fdofs, cdofs);
K_cf = K(cdofs, fdofs); K_cc = K(cdofs, cdofs);
a_f = a(fdofs); a_c = a(cdofs);
f_f = f(fdofs); f_c = f(cdofs);
```
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
Focusing on the matrix ``K_{ij}``, we note that if we split our domain into ``N_e`` elements, with start and endpoints ``a_i`` and ``b_i`` (where ``i`` is the element number), we can equally well evaluate this integral as a sum of the contributions for each element:
```math
K_{ij} = \int_{a}^{b} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x = \sum_{i = 1}^{N_e} \int_{a_i}^{b_i} \frac{\partial N_i}{\partial x}\ k\ \frac{\partial N_j}{\partial x}\ \mathrm{d}x
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
N_2^{(1)}(x) & x \in e_1 \\ 
N_1^{(2)}(x) & x \in e_2 
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
Even in the global numbering, we always have that the shape function corresponding to a specific node is 1 at that node, and all other shape functions are zero.Therefore, if we have Neumann boundary conditions, (known flux), at the right boundary, we want to add the contribution to the global load vector, specifically the contribution
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
f(right_node) = f(right_node) + qb
```
where `qb` is our known boundary flux. 

#### Dirichlet BC
Utilizing the same property, that ``N_i(x_j)`` is 1 if ``i = j``, and 0 if ``i \neq j``, we have that the temperature on the left side, ``T(a) = a_1``, such that we can enforce the dirichlet boundary conditions and solve the equation system via partitioning as 
```
cdofs = node_sets{"left"}
ac = [12] % Set the temperature on the left side to 12 degrees
fdofs = setdiff(1:ndofs, cdofs) % All other dofs
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
### Divergence and divergence theorem
Consider a vector field, 
```math
\underline{q}(\underline{x}) = \begin{bmatrix} q_1(\underline{x}) \\ q_2(\underline{x}) \end{bmatrix}, \quad \underline{x} = \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
```
around a body with volume ``V``:
"""

# ╔═╡ 58c1a7d8-03fb-409d-8589-46d2fa597f82
md" **TODO:** Add figure showing rectangle inside a complete body, ``\Omega``, to derive: ``\mathrm{div}(\underline{q}) = \partial q_i / \partial x_i``"

# ╔═╡ 8786dc4e-f7f4-4ccc-85f5-5302988c696c
md"""
The definition of the divergence of ``\underline{q}`` is
```math
\mathrm{div}(\underline{q}) := \lim_{V\rightarrow 0}\left[ \frac{1}{V} \int_{\Gamma} \underline{q} \cdot \underline{n}\ \mathrm{d}\Gamma \right]
```
**Some derivations...**
"""

# ╔═╡ 5608fe81-9b47-4214-a8ad-0960ca147b27
md"""
**TODO:** Use the same ``\Omega``, split into ``\Omega_1`` and ``\Omega_2``, to show that
```math
\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \int_{\Gamma_1} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma + \int_{\Gamma_2} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma
```
"""

# ╔═╡ 33618c32-a221-487c-b1ac-0b4beb347c7d
md"""
Now we can show that we can do
```math
\int_\Gamma \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma = \sum_{i = 1}^{N_V} \left[V_i \left[\frac{1}{V_i} \int_{\Gamma_i} \underline{q}\cdot\underline{n}\ \mathrm{d}\Gamma\right]\right]
```
And dividing the body into infinitely many subvolumes, i.e. letting ``V_i \rightarrow 0``, we get the divergence theorem
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

# ╔═╡ fccdcef3-3768-47d4-ba0e-ecca7308d356
md"""
### Weak form

"""

# ╔═╡ 27087e5b-08fb-4e0e-962c-30999860e238
md"""
### FE form

"""

# ╔═╡ d9385605-737f-401b-aba9-5d8d7a18a18c
md"""
## L7-8: Finite elements in 2D and 3D
### Shape functions
### Parametric elements
### Numerical integration
### Mapping of gradients
### Dirichlet boundary conditions
### Neumann boundary conditions
"""

# ╔═╡ 55570b97-a717-4a7c-89c4-a2ae941a8e32
md"""
## L9: Transient heat flow
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
## L12a: Mechanical equilibrium
"""

# ╔═╡ f13237a3-1bf4-4ec9-a797-20458df04e52
md"""
## L12b: Linear elasticity
"""

# ╔═╡ 5f0cf5a1-bd8c-4637-a28d-e4fd134254ab
md"""
## L13: Weak form of the mechanical equilibrium
"""

# ╔═╡ 976fbb59-1adf-4e95-84b8-75c7bf302917
md"""
## L14: Finite element analysis of linear elasticity
"""

# ╔═╡ f61e8d2d-1842-4ba8-aff8-01f6b9ad9025
md"""
## L15: Postprocessing and Robin BC
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Ferrite = "c061ca5d-56c9-439f-9c0e-210fe06d3992"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
CairoMakie = "~0.15.6"
Ferrite = "~1.1.0"
LaTeXStrings = "~1.4.0"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
manifest_format = "2.0"
project_hash = "3e77b724f1f6ec76c265ca194d0f1c0090a63c93"

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
git-tree-sha1 = "f8caabc5a1c1fb88bcbf9bc4078e5656a477afd0"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.6"

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
git-tree-sha1 = "cb1299fee09da21e65ec88c1ff3a259f8d0b5802"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.4"

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
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

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
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

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
version = "1.6.0"

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
git-tree-sha1 = "eaa040768ea663ca695d442be1bc97edfe6824f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.3+0"

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
git-tree-sha1 = "569ab58080263a9e8f7870e093817e3d8edafcc7"
uuid = "c061ca5d-56c9-439f-9c0e-210fe06d3992"
version = "1.1.0"

    [deps.Ferrite.extensions]
    FerriteBlockArrays = "BlockArrays"
    FerriteMetis = "Metis"

    [deps.Ferrite.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b"

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
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

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
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"
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
git-tree-sha1 = "ba6ce081425d0afb2bedd00d9884464f764a9225"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.2.2"
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
git-tree-sha1 = "bf0210c01fb7d67c31fed97d7c1d1716b98ea689"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.1"

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
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"

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
version = "8.11.1+1"

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
git-tree-sha1 = "368542cde25d381e44d84c3c4209764f05f4ef19"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.6"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "a370fef694c109e1950836176ed0d5eabbb65479"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.6"

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
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

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
version = "3.5.1+0"

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
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

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
git-tree-sha1 = "1f7f9bbd5f7a2e5a9f7d96e51c9754454ea7f60b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.4+0"

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
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

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
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

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
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

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
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

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
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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
version = "5.13.1+1"

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
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"

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
# ╟─33618c32-a221-487c-b1ac-0b4beb347c7d
# ╟─4b7acd47-9515-4325-bfb3-2899f847c26b
# ╟─b1999d4a-3111-4203-ab1c-173905f06d80
# ╟─c1e878bd-0aac-4f28-a528-75c9ce9e3e4a
# ╠═fccdcef3-3768-47d4-ba0e-ecca7308d356
# ╠═27087e5b-08fb-4e0e-962c-30999860e238
# ╟─d9385605-737f-401b-aba9-5d8d7a18a18c
# ╠═55570b97-a717-4a7c-89c4-a2ae941a8e32
# ╟─509c3307-b814-4da9-aecc-7adb47f8ec95
# ╠═1f5672bc-28af-4ab7-b159-197ebf1e12a3
# ╠═f13237a3-1bf4-4ec9-a797-20458df04e52
# ╠═5f0cf5a1-bd8c-4637-a28d-e4fd134254ab
# ╠═976fbb59-1adf-4e95-84b8-75c7bf302917
# ╠═f61e8d2d-1842-4ba8-aff8-01f6b9ad9025
# ╟─7e7cf116-f844-4345-bb00-e9829e6262be
# ╟─6bfc5b5e-cb6f-4e33-8b86-f99ffa713eb5
# ╟─f436cc3b-49cb-4f7b-915e-e2c0fee4c2cb
# ╟─87f373ca-78eb-4e1f-b539-33d3fd05ab21
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
