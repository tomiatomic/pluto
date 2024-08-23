### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 07e9d298-bc12-4007-ad79-43dd2e4dccdd
using PlutoUI

# ╔═╡ 3ed350e4-992e-42ae-bca3-04229e8cd364
TableOfContents()

# ╔═╡ 02a6191a-cd7b-402c-b3e1-7db1c6ebb36b
md"# Quasiparticle spectra of *2H-NbSe₂*: Two-band superconductivity and the role of tunneling selectivity
[Noat et al. PRB2015](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.134510)"

# ╔═╡ 80cbc3b4-6bc7-4d0a-a769-92b5c5a99a14
md"## equation (1)
$\Delta_i(E) = \frac{\Delta_i^0 + \Gamma_{ij}\Delta_j(E)/\sqrt{\Delta_j^2(E)-E^2}}
{1+\Gamma_{ij}\sqrt{\Delta_j^2(E)-E^2}}$"

# ╔═╡ 6fa5d977-f472-4454-9e58-4578e0810ea5
md"## equation (2)
$N_S^i(E) = N_i(E_F)\Re\left[\frac{|E|}
{\sqrt{E^2-\Delta_i(E)^2}}\right]$"

# ╔═╡ 52a4dd88-7c7a-4d67-b51a-be2901275187
md"# Thèse de Ronan Lamy
Here's the link to [theses.fr](https://theses.fr/2004PA066555)"

# ╔═╡ ef41f6ab-fbe8-4f8f-a2fb-4ee310d09d7e
md"## modèle de McMillan
Here's the link to [McMillan PR1968](https://journals.aps.org/pr/abstract/10.1103/PhysRev.175.537)"

# ╔═╡ bc969468-7497-4521-b03f-7830aa9faeb4
md"$\Delta_\kappa(i\omega_n) = 
\left
(\Delta_\kappa^0 + 
\sum\limits_{\lambda \neq \kappa}\Gamma_{\kappa,\lambda}{ \frac{\Delta_\lambda(i\omega_n)}
{\sqrt{\omega_n^2 + \Delta_\lambda^2(i\omega_n)}}}
\right) \bigg/
\left
(1 + 
\sum\limits_{\lambda \neq \kappa}\Gamma_{\kappa,\lambda}{ \frac{1}
{\sqrt{\omega_n^2 + \Delta_\lambda^2(i\omega_n)}}}
\right) ~(1.71)$
*remplacer $\sqrt{\omega_n^2 + \Delta_\lambda^2(i\omega_n)}$ par $-i\sqrt{\omega_n^2 - \Delta_\lambda^2(\omega)}...$"

# ╔═╡ 579c11be-0386-4941-9e36-55d6014bcc50
md"## équation differentielle
$\begin{align}
&\Delta_1'(\omega) = \\
&\Bigg[\omega(\Delta_1-\Delta_1^0)\Bigg((\Delta_1-\Delta_2)(\Delta_2-\Delta_2^0)^3\Gamma_{1,2}^2+(\Delta_1-\Delta_1^0)^2\bigg(\Delta_2(\Delta_2-\Delta_2^0)^3 \\
&+(-\Delta_1+\Delta_2)(\Delta_1-\Delta_2^0)\Gamma_{2,1}^2\bigg)\Bigg)\Bigg]\Bigg/ \\
&\Bigg/\Bigg[\Delta_2(\Delta_1-\Delta_1^0)^3(\Delta_2-\Delta_2^0)\bigg(\Delta_1(\Delta_2-\Delta_2^0)^2-(\Delta_1-\Delta_2)\Gamma_{2,1}^2\bigg) \\
&+(\Delta_1-\Delta_2)\Gamma_{1,2}^2\bigg(\Delta_1(\Delta_1-\Delta_1^0)(\Delta_2-\Delta_2^0)^3-(\Delta_1-\Delta_2)^2(\Delta_1^0-\Delta_2^0)\Gamma_{2,1}^2\bigg)\Bigg]
\end{align}~(1.72)$

#### or

$\begin{gather}
\Delta'_i(\omega)= \Bigl[\omega(\Delta_i-\Delta_i^0)\\
\left((\Delta_i-\Delta_j)(\Delta_j-\Delta_j^0)^3\Gamma_{1,2}^2+(\Delta_i-\Delta_i^0)^2
\left(\Delta_j(\Delta_j-\Delta_j^0)^3+(-\Delta_i+\Delta_j)(\Delta_i-\Delta_j^0)\Gamma_{j,i}^2\right)\right)\Bigr]\\
\bigg/\Bigl[\Delta_j(\Delta_i-\Delta_i^0)^3(\Delta_j-\Delta_j^0)\left(\Delta_i(\Delta_j-\Delta_j^0)^2-(\Delta_i-\Delta_j)\Gamma_{j,i}^2 \right)\\
+(\Delta_i-\Delta_j)\Gamma_{i,j}^2\left(\Delta_i(\Delta_i-\Delta_i^0)(\Delta_j-\Delta_j^0)^3-(\Delta_i-\Delta_j)^2(\Delta_i^0-\Delta_j^0)\Gamma_{j,i}^2\right)
\Bigr]
\end{gather}$
"

# ╔═╡ f14ad140-67df-4fc1-8b57-39d0ba6d84da
md"""# $\LaTeX$ of "McMillan 2020.nb" 
Mathematica functions trnascribed to Julia"""

# ╔═╡ e3f0afb0-c5dd-11ee-1a99-a53376ffbbcc
md"## recursiondelta (500 $\times$ by deltaini)
$\begin{gather}
x = \frac{{\Delta_1 \sqrt{{- (y^2 - E^2)}} + i \Gamma_1 y}}{{\sqrt{{- (y^2 - E^2)}} + i \Gamma_1}} \\
\\
y = \frac{{\Delta_2 \sqrt{{- (x^2 - E^2)}} + i \Gamma_2 x}}{{\sqrt{{- (x^2 - E^2)}} + i \Gamma_2}}
\end{gather}$"

# ╔═╡ 9c189f39-4c17-4ea8-85c6-352b1a2a4abe
md"## deltadif
$\begin{align}
\Delta_1^2 u &- \left(\Delta_1 u^2 + y[u] (\Gamma_1^2 - \Delta_1 y[u])\right) x'[u] \\
&+ (-\Delta_1^2 + \Gamma_1^2) y[u] y'[u] + x[u]^2 (u - y[u]y'[u])\\ 
&+ x[u] \left(-2 \Delta_1 u + (\Gamma_1^2 + u^2 - y[u]^2) x'[u] - (\Gamma_1^2 - 2 \Delta_1 y[u])y'[u]\right) = 0
\end{align} ~~~(1)$
$\begin{align}
\Delta_2^2 u &- \left(\Delta_2 u^2 + x[u] (\Gamma_2^2 - \Delta_2 x[u])\right) y'[u] \\
&+ (-\Delta_2^2 + \Gamma_2^2) x[u] x'[u] + y[u]^2 (u - x[u] x'[u]) \\
&+ y[u] \left(-2 \Delta_2 u + (\Gamma_2^2 + u^2 - x[u]^2) y'[u]  - (\Gamma_2^2 - 2 \Delta_2 x[u]) x'[u]\right) = 0
\end{align} ~~~(2)$
#### or
$\begin{align}
\Delta_1^2 u &- \left(\Delta_1 u^2 + y[u] (\Gamma_1^2 - \Delta_1 y[u])\right) \frac{dx}{du} \\
&+ (-\Delta_1^2 + \Gamma_1^2) y[u] \frac{dy}{du} + x[u]^2 (u - y[u] \frac{dy}{du}) \\
&+ x[u] \left(-2 \Delta_1 u + (\Gamma_1^2 + u^2 - y[u]^2) \frac{dx}{du} - (\Gamma_1^2 - 2 \Delta_1 y[u]) \frac{dy}{du}\right) = 0
\end{align} ~~~(1)$
$\begin{align}
\Delta_2^2 u &- \left(\Delta_2 u^2 + x[u] (\Gamma_2^2 - \Delta_2 x[u])\right) \frac{dy}{du} \\
&+ (-\Delta_2^2 + \Gamma_2^2) x[u] \frac{dx}{du} + y[u]^2 (u - x[u] \frac{dx}{du}) \\
&+ y[u] \left(-2 \Delta_2 u + (\Gamma_2^2 + u^2 - x[u]^2) \frac{dy}{du} - (\Gamma_2^2 - 2 \Delta_2 x[u]) \frac{dx}{du}\right)  = 0
\end{align} ~~~(2)$
"

# ╔═╡ d939159c-19b6-4a46-aa0c-ba4e9ef19b1a
md"## nMcMillan
$N_i[E] = \Re\frac{E}{\sqrt{E^2-\Delta_iE^2}}$"

# ╔═╡ e39e7261-ac94-480c-aea0-b7216e1cd469
md"# My Pluto notebook"

# ╔═╡ b4c5c459-8d1d-43d1-a47e-ebf49a48deb9
md"## recdel (500 $\times$ by delini)
$\begin{gather}
x = \frac{{\Delta_1 \sqrt{{- (y^2 - E^2)}} + i \Gamma_1 y}}{{\sqrt{{- (y^2 - E^2)}} + i \Gamma_1}} \\
\\
y = \frac{{\Delta_2 \sqrt{{- (x^2 - E^2)}} + i \Gamma_2 x}}{{\sqrt{{- (x^2 - E^2)}} + i \Gamma_2}}
\end{gather}$"

# ╔═╡ 6e71bd4e-3736-4a9c-a064-7741a7bb70f3
md"## deldif
$\begin{gather}
\frac{dx}{dt} = \frac{\frac{dy}{dt}\left((-\Delta_1^2 + \Gamma_1^2)y - x^2y - (\Gamma_1^2 - 2\Delta_1y)x\right) + t(\Delta_1^2 - 2\Delta_1x + x^2)}{\left(\Delta_1t^2 + y(\Gamma_1^2 - \Delta_1y)\right) - (\Gamma_1^2 + t^2 - y^2)x}\\
\\
\frac{dy}{dt} = \frac{\frac{dx}{dt}\left((-\Delta_2^2 + \Gamma_2^2)x - y^2x - (\Gamma_2^2 - 2\Delta_2x)y\right) + t(\Delta_2^2 - 2\Delta_2y + y^2)}{\left(\Delta_2t^2 + x(\Gamma_2^2 - \Delta_2x)\right) - (\Gamma_2^2 + t^2 - x^2)y}
\end{gather}$"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.55"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "f64cdffc70331b0a2f407efefd54fd84eb680773"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─07e9d298-bc12-4007-ad79-43dd2e4dccdd
# ╟─3ed350e4-992e-42ae-bca3-04229e8cd364
# ╟─02a6191a-cd7b-402c-b3e1-7db1c6ebb36b
# ╟─80cbc3b4-6bc7-4d0a-a769-92b5c5a99a14
# ╟─6fa5d977-f472-4454-9e58-4578e0810ea5
# ╟─52a4dd88-7c7a-4d67-b51a-be2901275187
# ╟─ef41f6ab-fbe8-4f8f-a2fb-4ee310d09d7e
# ╟─bc969468-7497-4521-b03f-7830aa9faeb4
# ╟─579c11be-0386-4941-9e36-55d6014bcc50
# ╟─f14ad140-67df-4fc1-8b57-39d0ba6d84da
# ╟─e3f0afb0-c5dd-11ee-1a99-a53376ffbbcc
# ╟─9c189f39-4c17-4ea8-85c6-352b1a2a4abe
# ╟─d939159c-19b6-4a46-aa0c-ba4e9ef19b1a
# ╟─e39e7261-ac94-480c-aea0-b7216e1cd469
# ╟─b4c5c459-8d1d-43d1-a47e-ebf49a48deb9
# ╟─6e71bd4e-3736-4a9c-a064-7741a7bb70f3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
