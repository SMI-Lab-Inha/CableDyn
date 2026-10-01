.. SPDX-License-Identifier: Apache-2.0

Theory
======

This chapter states the equations CableDyn solves: the kinematics and the two line
elements, the constitutive laws, the hydrodynamic and seabed loads, the static and dynamic
solution procedures, and the models for attached objects. Units, frames, sign rules, and the
tension definition are fixed in :doc:`conventions` and are not repeated here. The
element-level implementation of each path, and the secondary Cosserat path, are described in
:doc:`solver_paths`.

Symbols are defined in the section that uses them. A few letters are reused: :math:`\beta` is the
wave direction, the Syrope fast-stiffness slope, a Newmark parameter, and an eigenvalue
parameter of the sagging-cable modes; :math:`T` is a tension and, in the wave formulas, a
period.

.. contents:: On this page
   :local:
   :depth: 2

Kinematics and discretisation
-----------------------------

A line is a curve :math:`\mathbf{r}(s, t) \in \mathbb{R}^3`, parametrised by the unstretched
arc length :math:`s \in [0, L]`. It is built from sections joined end to end; each section is
meshed into two-node elements of unstretched length :math:`L_0`, with the local coordinate
:math:`\xi = (s - s_e)/L_0 \in [0, 1]` on element :math:`e`. The axial strain and the
centreline curvature are the exact large-deformation measures

.. math::

   \varepsilon = \lvert \mathbf{r}' \rvert - 1, \qquad
   \kappa = \frac{\lvert \mathbf{r}' \times \mathbf{r}'' \rvert}{\lvert \mathbf{r}' \rvert^{3}},
   \qquad (\cdot)' = \partial(\cdot)/\partial s .

Neither measure is linearised, so both elements are geometrically exact for arbitrary
displacement and rotation. Which element a section uses is set by the bending stiffness of
its line type: ``EI = 0`` selects the tension element, ``EI > 0`` the bending element.

The EI = 0 tension element
~~~~~~~~~~~~~~~~~~~~~~~~~~

Each node carries its position only (three degrees of freedom). For nodal positions
:math:`\mathbf{r}_a, \mathbf{r}_b`,

.. math::

   \boldsymbol{\delta} = \mathbf{r}_b - \mathbf{r}_a, \quad
   \ell = \lvert\boldsymbol{\delta}\rvert, \quad
   \mathbf{t} = \boldsymbol{\delta}/\ell, \quad
   \varepsilon = \ell/L_0 - 1, \quad
   T = EA\,\varepsilon .

The internal force and its exact tangent (material plus geometric stiffness) are

.. math::

   \mathbf{f}_\text{int} = \begin{bmatrix} -T\,\mathbf{t} \\ +T\,\mathbf{t} \end{bmatrix},
   \qquad
   \mathbf{K}_t = \begin{bmatrix} \mathbf{B} & -\mathbf{B} \\ -\mathbf{B} & \mathbf{B} \end{bmatrix},
   \qquad
   \mathbf{B} = \frac{EA}{L_0}\,\mathbf{t}\mathbf{t}^{\mathsf T}
              + \frac{T}{\ell}\left(\mathbf{I} - \mathbf{t}\mathbf{t}^{\mathsf T}\right).

In **tension-only** mode an element with :math:`\varepsilon < 0` carries zero elastic force and
zero material stiffness, so force and tangent remain consistent. The deck and coupled routes
march the dynamics tension-only: a slack chain or rope exerts no elastic force. Axial damping
(``BA``, below) stays active on a slack element, as in MoorDyn, so a slack element with
:math:`BA \ne 0` can carry a damping force while it stretches or shortens. The static initial
condition is the tension-only equilibrium of the same law (see
:ref:`static-ei0`). The consistent mass of an element is

.. math::

   \mathbf{M}_e = \frac{\rho_A L_0}{6}
   \begin{bmatrix} 2\mathbf{I} & \mathbf{I} \\ \mathbf{I} & 2\mathbf{I} \end{bmatrix},

with :math:`\rho_A` the dry mass per unstretched metre.

The cubic-Hermite bending element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each node carries its position and a material tangent, :math:`[\mathbf{r},\ \mathbf{m}]`
with :math:`\mathbf{m} = \mathbf{r}'`, so an element has 12 degrees of freedom and no rotation
variables. The centreline is the cubic-Hermite interpolant

.. math::

   \mathbf{r}(\xi) = H_1(\xi)\,\mathbf{r}_1 + H_2(\xi)\,\mathbf{m}_1
                   + H_3(\xi)\,\mathbf{r}_2 + H_4(\xi)\,\mathbf{m}_2,

   H_1 = 1 - 3\xi^2 + 2\xi^3,\quad H_2 = L_0(\xi - 2\xi^2 + \xi^3),\quad
   H_3 = 3\xi^2 - 2\xi^3,\quad H_4 = L_0(\xi^3 - \xi^2),

which is :math:`C^1` across nodes, so curvature is defined pointwise inside every element.
The magnitude :math:`\lvert\mathbf{m}\rvert = 1 + \varepsilon` carries the axial stretch at a
node; only its direction is a bending quantity. The element stores the axial-plus-bending
strain energy

.. math::

   U_e = \int_0^{L_0} \left( \tfrac{1}{2} EA\,\varepsilon^2 + \tfrac{1}{2} EI\,\kappa^2 \right)
         \mathrm{d}s ,
   \qquad
   \kappa^2 = \frac{\lvert \mathbf{a} \times \mathbf{b} \rvert^2}{\lvert\mathbf{a}\rvert^6},
   \quad \mathbf{a} = \mathbf{r}',\ \mathbf{b} = \mathbf{r}'' .

Using :math:`\kappa^2` directly, rather than :math:`\kappa`, keeps the energy smooth through
the straight state. The stress-free reference is straight, and the reported bending moment is
:math:`M = EI\,\kappa`.

At each Gauss point, :math:`\mathbf{a}` and :math:`\mathbf{b}` are fixed linear maps of the
12 element DOFs. The gradient and Hessian of the energy density with respect to
:math:`(\mathbf{a}, \mathbf{b})` are written in closed form and chained to the DOFs, giving
:math:`\mathbf{f}_\text{int} = \partial U/\partial \mathbf{q}` and the symmetric tangent
:math:`\mathbf{K}_t = \partial^2 U/\partial \mathbf{q}^2` exactly. Automatic differentiation of
the same energy is used only in the test suite, as an independent check of the closed form,
together with finite differences of :math:`\mathbf{f}_\text{int}`.

**Quadrature.** Axial and bending energies use Gauss–Legendre integration, four points each by
default. The two orders are independently selectable from 1 to 6 (options
``axial_quadrature_order`` and ``bending_quadrature_order``); unequal orders give a selective
integration. The consistent mass
:math:`\mathbf{M}_e = \rho_A \int_0^{L_0} \mathbf{N}^{\mathsf T}\mathbf{N}\,\mathrm{d}s` uses the same
Hermite shapes and a four-point rule, which is exact for the degree-6 integrand; it is
block-diagonal per Cartesian component. Distributed hydrodynamic loads use three-point Gauss
integration.

**Curvature recovery.** Curvature outputs are evaluated from the interpolated centreline, not
only at nodes or Gauss points. The peak curvature of an element is located by sampling 16
equal intervals and bisecting each at least five and at most twelve times, continuing where the
curvature departs from linear variation by more than 5 % of the local maximum. The same search
locates the extreme axial resultant :math:`N = EA(\lvert\mathbf{r}'\rvert - 1)` inside every
element, with a 2 % criterion relative to the largest local magnitude (at least 1 N).

The element follows the geometrically exact Kirchhoff-rod formulation of
:ref:`Boyer et al. (2011) <ref-boyer2011>` with the :math:`C^1` Hermite centreline
interpolation of :ref:`Meier, Popp & Wall (2015) <ref-meier2015>`, reduced to a torsion-free
cable. Torsion enters only through the condensed line term below.

.. _theory-torsion:

Condensed torsion
~~~~~~~~~~~~~~~~~

A line restrained in torsion at both ends (``END CONNECTIONS`` ``TorsStiffness``, see
:doc:`driver_format`) adds the twist of an isotropic Kirchhoff rod without adding degrees of
freedom. For an isotropic section with no distributed torque the twisting moment
:math:`M = GJ\,u_3`, with :math:`u_3` the material twist rate, is uniform along the rod
(:ref:`van der Heijden et al., 2003 <ref-vanderheijden2003>`). The twist field
then condenses to one scalar per line, and the line energy gains

.. math::

   E_t = \frac{(\Phi - \Theta(\mathbf{q}))^2}{2C}, \qquad
   M = \frac{\Phi - \Theta}{C}, \qquad
   C = \sum_e \frac{L_{0,e}}{GJ_e} + \frac{1}{k_A} + \frac{1}{k_B} .

:math:`\Phi` is the imposed relative roll of the End B frame with respect to the End A frame
(``Pretwist`` and the ``motionFile`` roll column), :math:`k_A, k_B` are the optional torsional end
springs (:math:`1/k = 0` for ``Rigid``), and :math:`\Theta(\mathbf{q})` is the geometric twist of
the centreline: the angle about the End B direction :math:`\mathbf{d}_B`, right-handed, from the
End B reference normal :math:`\mathbf{n}_B` to the End A normal :math:`\mathbf{n}_A` carried to End B
by parallel transport. Each end frame :math:`(\mathbf{d}, \mathbf{n})` turns with the body or vessel
that holds the end. A clamped end has :math:`\mathbf{d}` equal to its end tangent; at a
bending-pinned end :math:`\mathbf{d}` is the connection direction and the transport starts with
the smallest rotation from :math:`\mathbf{d}` to the tangent.

Because the Hermite centreline is :math:`C^1`, the transport factorises into one holonomy per
element and a chain of smallest rotations between the unit node tangents
:math:`\mathbf{t}_k = \mathbf{m}_k/\lvert\mathbf{m}_k\rvert`:

.. math::

   \Theta = \sum_e h_e + \chi, \qquad
   h_e = \int_e \frac{\mathbf{t}_1 \cdot (\mathbf{a} \times \mathbf{b})}
                     {\lvert\mathbf{a}\rvert^2 \left(1 + \mathbf{t}_1\cdot\mathbf{a}/\lvert\mathbf{a}\rvert\right)}
         \,\mathrm{d}s ,

with :math:`\mathbf{t}_1` the tangent at the element's first node, :math:`\mathbf{a} = \mathbf{r}'`,
:math:`\mathbf{b} = \mathbf{r}''` at the Gauss points of the bending rule, and :math:`\chi` the
angle of the smallest-rotation chain
:math:`\mathbf{d}_A \to \mathbf{t}_1 \to \dots \to \mathbf{t}_N \to \mathbf{d}_B` applied to
:math:`\mathbf{n}_A`. For differentiation, each link of the chain is written with a gauge frozen at
the current configuration as a sum of signed solid angles
:math:`\Omega(\mathbf{x}, \mathbf{y}, \mathbf{z}) = 2\operatorname{atan2}(\mathbf{x}\cdot(\mathbf{y}\times\mathbf{z}),
1 + \mathbf{x}\cdot\mathbf{y} + \mathbf{y}\cdot\mathbf{z} + \mathbf{z}\cdot\mathbf{x})` of geodesic
triangles, each of which involves at most two neighbouring tangents. The gradient and Hessian
of :math:`\Theta` are therefore local, in closed form, and have the half-bandwidth of the cable
tangent; the value of :math:`\Theta` is taken from the sequential product. This is the
discrete parallel transport and holonomy of :ref:`Bergou et al. (2008) <ref-bergou2008>` applied
to the Hermite centreline. The test suite checks the closed form against automatic
differentiation of an independent reference implementation, against finite differences, and
for invariance under rigid rotation.

:math:`\Theta` is known modulo :math:`2\pi`. The solver carries its unwrapped value as a state
of the line, unwraps each new evaluation to the branch nearest the last accepted value, and
accepts no step that changes :math:`\Theta` by more than :math:`\pi/2`; such a step is cut. A
link or Gauss point at which the tangent turns by more than 120° within one element or between
neighbouring nodes (:math:`1 + \mathbf{t}_1\cdot\mathbf{t} < 0.5`) stops with a named error; a
mesh that resolves the curvature is far from it.

**Residual and tangent.** The torsion term adds :math:`-M\,\nabla\Theta` to the residual and

.. math::

   \mathbf{K} = \underbrace{\mathbf{K}_c - M\,\nabla^2\Theta}_{\mathbf{B}\ \text{(banded)}}
               + \frac{1}{C}\,\nabla\Theta\,\nabla\Theta^{\mathsf T}

to the tangent. The rank-one term is not banded; each Newton step solves the bordered system
by Sherman–Morrison on the band factorisation of :math:`\mathbf{B}` with two right-hand sides,
:math:`\mathbf{B}\mathbf{y} = -\mathbf{R}`, :math:`\mathbf{B}\mathbf{z} = \nabla\Theta`,
:math:`\mathbf{x} = \mathbf{y} - \mathbf{z}\,(\nabla\Theta\cdot\mathbf{y})/(C + \nabla\Theta\cdot\mathbf{z})`,
and accepts it on a backward-error test, with a shifted factorisation and iterative refinement
as the fallback. Dropping the rank-one term would leave a Newton iteration that converges only
linearly, at the rate :math:`\nabla\Theta^{\mathsf T}\mathbf{B}^{-1}\nabla\Theta / C`; on the
post-buckled states of the validation (see :doc:`validation`) that rate is 0.05 to 0.22. The
same term enters the dynamic effective tangent, scaled by :math:`1 - \alpha_f` of the force
blend, and is frozen together with a reused factorisation.

**Statics and stability.** The imposed twist is the last static load stage, ramped from the
geometric twist of the untwisted equilibrium in steps of at most :math:`\pi/4`. Each converged
stage is tested for stability by the inertia of :math:`\mathbf{K}` on the free DOFs,
:math:`\operatorname{neg}(\mathbf{K}) = \operatorname{neg}(\mathbf{B}) -
[1 + \nabla\Theta^{\mathsf T}\mathbf{B}^{-1}\nabla\Theta/C < 0]`, with the negative eigenvalues
of :math:`\mathbf{B}` counted by a banded :math:`\mathbf{L}\mathbf{D}\mathbf{L}^{\mathsf T}`
factorisation. Above a buckling onset the straight twisted state is a saddle to which Newton
converges quadratically; the static solve then descends along the lowest mode with energy
acceptance onto the buckled branch, and reports the descent. The onset is that of a
clamped–clamped rod under the dead tension :math:`T`, eq. (33) of
:ref:`van der Heijden et al. (2003) <ref-vanderheijden2003>`; at :math:`T = 0` it reduces to
:math:`\tan x = x`, :math:`x = ML/2EI`, so :math:`M_{cr} = 8.9868\,EI/L`. At large tension it
approaches the localised value :math:`2\sqrt{EI\,T}` of an infinite rod from above, reported in
the same paper. A bending-pinned end restrained in torsion transmits the torque about the
bisector of the connection direction and the tangent (a *semi-tangential* end, a
constant-velocity joint), for which :math:`\tan x = -x/3` and :math:`M_{cr} = 4.9113\,EI/L` at
zero tension, not the :math:`2\pi EI/L` of Greenhill's axial-torque hinge. OrcaFlex's
bending-free end with a twisting spring gives the same value (see :doc:`validation`).

**Dynamics and loads.** The torsion force enters the force-blended residual like the internal
force, so a dynamic deck with torsion requires ``True alpha_force_blend`` (the default). Torsion
carries no inertia: the torque follows the imposed twist within the step, without a torsional
wave. This is accurate while the first torsional frequency of the line,
:math:`\sqrt{GJ/I_p}/(2L)` with :math:`I_p` the polar mass moment per length, is well above the
excitation; it is not for very long or torsionally soft lines. The moment of the torque on the
object at End A is :math:`M\,\partial\Theta/\partial\boldsymbol{\omega}_A`, the derivative with
respect to a rotation :math:`\boldsymbol{\omega}_A` of that object (for a clamped end, of its frame
and the end tangent together), so the virtual work of the body moment is exact; for a straight
line it is the torque along the line axis. The static body equilibrium uses the same moment and
its stiffness :math:`\mathbf{a}\mathbf{a}^{\mathsf T}/C`, :math:`\mathbf{a} = \partial\Theta/\partial\boldsymbol{\omega}_A`.

**Limits.** Seabed friction does not resist twist, so the laid part of a line twists freely
(the friction torque of a cable on the seabed can be comparable with the twist torque near
touchdown). There is no torque–tension coupling of armoured cables or wire ropes, no
anisotropic or pre-twisted section, and no self-contact: a loop that closes on itself is not
resolved. The model is that of :ref:`Meier, Popp & Wall (2014) <ref-meier2014>` and
:ref:`Bergou et al. (2008) <ref-bergou2008>` restricted to the isotropic, torque-free-span case;
it is not a new formulation.

Constitutive laws
-----------------

Axial response is chosen per line type. Stateful laws advance their committed internal state
once per converged step, from the converged kinematics, in a single commit that is rolled back
together with the positions if any part of the step fails.

Linear elasticity
~~~~~~~~~~~~~~~~~

:math:`T = EA\,\varepsilon`. On the ``EI = 0`` path the dynamic law is tension-only, as above.
On the Hermite path the axial law is two-sided: the energy
:math:`\tfrac{1}{2}EA\,\varepsilon^2` is used for both signs of :math:`\varepsilon`, which keeps the
bending equations well posed in locally slack regions.

**Tensile audit (Hermite path).** Option ``tensile_safety`` inspects the axial force of the
line after each converged step. The audited force is the one that enters equilibrium, as in
the static branch audit below: the element-mean axial force (the three-point Gauss mean of the
signed resultant :math:`EA\,(|\mathbf r'| - 1)`), averaged over each element and its two
neighbours weighted by length (one neighbour at a line end). The pointwise resultant of a
stiff-``EA`` cubic-Hermite element oscillates about that mean, and a single element mean dips
next to a nodal seabed-contact reaction, by amounts that grow with the element length: judged
pointwise, a coarse mesh whose tensions and curvatures are converged reads as compressed (the
Gulf of Mexico 80 m lazy-wave cable of :ref:`Lozon et al. (2025) <ref-lozon2025>` under a
3 m, 12 s heave: tens of kN on 48 elements, where the fairlead
tension and the peak curvature are within 1 % of a 1024-element mesh). A compressed span of
three or more elements is never averaged away. A step is flagged when the three-element mean
falls below :math:`-\epsilon_\text{tol}\,\overline{EA}`, with :math:`\overline{EA}` the
length-weighted ``EA`` of the same elements and :math:`\epsilon_\text{tol} = 2\times10^{-6}` by
default (``tensile_strain_tolerance``). In ``error`` mode the step is rejected; in ``warn``
mode it is committed and the worst event is reported (force, threshold, centre element, time).
The audit is a response diagnostic and does not change the constitutive law.

Axial damping
~~~~~~~~~~~~~

Axial Kelvin–Voigt damping adds a resultant proportional to the strain rate. On the
``EI = 0`` path, for an element with tangent :math:`\mathbf{t}` and nodal velocities
:math:`\mathbf{v}_a, \mathbf{v}_b`,

.. math::

   T_d = \frac{BA}{L_0}\,\mathbf{t}\cdot(\mathbf{v}_b - \mathbf{v}_a),

applied along :math:`\mathbf{t}` with the MoorDyn sign convention
(:ref:`Hall & Goupee, 2015 <ref-hall2015>`). On the Hermite path the same law is integrated
over the element,

.. math::

   \dot\varepsilon = \frac{\mathbf{t}\cdot\mathbf{v}_{,\xi}}{L_0}, \qquad
   \mathbf{t} = \frac{\mathbf{r}_{,\xi}}{\lvert\mathbf{r}_{,\xi}\rvert}, \qquad
   \mathbf{f}_{d,a} = \int_0^1 BA\,\dot\varepsilon\,H'_a(\xi)\,\mathbf{t}\,\mathrm{d}\xi ,

with exact position and velocity Jacobians. A negative ``BA`` column value is a target
damping ratio :math:`\zeta = -BA_\text{input}` of the axial mode of each segment, resolved per
element as

.. math::

   BA = \zeta\,L_0\,\sqrt{EA\,\rho_A},

using the final element length after any mesh refinement. The damping share is included in
the reported tension.

Viscoelastic rope (MoorDyn ``ElasticMod`` 2 and 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viscoelastic law of :ref:`Hall, Duong & Lozon (2023) <ref-hall2023>` is a four-parameter
series-Kelvin solid: a fast Kelvin–Voigt branch (spring :math:`EA_D` parallel to dashpot
:math:`BA_D`) in series with a slow Kelvin–Voigt branch (spring :math:`EA_1` parallel to dashpot
:math:`BA_s`). It reduces to the standard linear solid only when :math:`BA_D = 0`. The slow
spring is chosen so that the static composite stiffness equals the line-type ``EA``:

.. math::

   EA_1 = \frac{EA_D\,EA}{EA_D - EA}, \qquad EA_D > EA .

Each element carries the slow-branch stretch :math:`\Delta l_1` as a state. With
:math:`\Delta l = \ell - L_0`,

.. math::

   \dot{\Delta l_1} = \frac{EA_D\,\Delta l - (EA_D + EA_1)\,\Delta l_1 + BA_D\,\dot\ell}
                           {BA_D + BA_s},
   \qquad
   T = \frac{EA_1\,\Delta l_1 + BA_s\,\dot{\Delta l_1}}{L_0} .

The series chain is tension-only: the element carries :math:`\max(T, 0)`, and while
:math:`T = 0` the slow branch relaxes under zero load, :math:`BA_s\,\dot{\Delta l_1} + EA_1\,\Delta l_1 = 0`.
At :math:`T = 0` the loaded and slack states coincide, so the force and the state are continuous
through slack. MoorDyn's formulation applies the slack condition to the elastic part when
:math:`\Delta l < 0`, so the two formulations differ through slack while
:math:`\Delta l_1 \ne 0` and agree whenever MoorDyn's tension is non-negative with
:math:`\Delta l \ge 0`. Inside a time step
the state equation is linear in :math:`\Delta l_1`, so backward Euler eliminates it in closed form,

.. math::

   \Delta l_1(\mathbf{q}, \mathbf{v}) =
   \frac{\Delta l_{1,n} + \Delta t\,(EA_D\,\Delta l + BA_D\,\dot\ell)/B}
        {1 + \Delta t\,(EA_D + EA_1)/B}, \qquad B = BA_D + BA_s ,

giving an element load with exact position and velocity Jacobians inside the implicit step.
The initial state is the zero-rate partition
:math:`\Delta l_1 = \Delta l\,(EA_D - EA)/EA_D`, at which the tension equals the static
:math:`EA\,\Delta l/L_0`, so the static equilibrium is an exact fixed point of the dynamics.

``ElasticMod`` 3 makes the dynamic stiffness depend on the mean load. For
:math:`\Delta l_1 > 0` and :math:`k = EA/L_0`,

.. math::

   EA_D = \tfrac{1}{2}\left(\alpha_\text{MBL} + v_\beta\,\Delta l_1\,k + EA
   + \sqrt{\alpha_\text{MBL}^2 + 2\alpha_\text{MBL}\,k\,(v_\beta\,\Delta l_1 - L_0)
   + k^2\,(v_\beta\,\Delta l_1 + L_0)^2}\right),

and :math:`EA_D = \alpha_\text{MBL}` for :math:`\Delta l_1 \le 0`. It is evaluated from the
committed state and held over the step. Parameter sets with :math:`EA_D \le EA`, a stiffness
condition :math:`(EA_D + EA)/(EA_D - EA) > 10^3`, negative dashpots, or :math:`BA_D + BA_s = 0`
fail closed. The viscoelastic law applies to ``EI = 0`` lines.

Syrope polyester model
~~~~~~~~~~~~~~~~~~~~~~

The Syrope model (:ref:`Falkenberg, Åhjem & Yang, 2017 <ref-falkenberg2017>`), in its
MoorDyn-F/C form, describes polyester from an **original working curve** (OWC), a
strain–tension table :math:`\varepsilon_\text{owc}(T)`. The strain is split into a fast part
with the load-dependent stiffness

.. math::

   \varepsilon_\text{fast}(T) = \frac{1}{\beta}\ln\!\left(1 + \frac{\beta}{\alpha}\,T\right),
   \qquad \frac{\mathrm{d}T}{\mathrm{d}\varepsilon_\text{fast}} = EA_d = \alpha + \beta\,T,

and a slow part :math:`\varepsilon_\text{slow} = \varepsilon_\text{curve}(T) - \varepsilon_\text{fast}(T)`,
where :math:`\varepsilon_\text{curve}` is the active curve (the OWC or the working curve below).
:math:`T_\text{curve}(\varepsilon_\text{slow})` denotes the inverse of this slow-strain relation:
the mean tension at which the active curve has the slow strain :math:`\varepsilon_\text{slow}`.
Loading beyond the running maximum mean tension :math:`T_\text{max}` follows the OWC. Below it
the line follows a **working curve** regenerated from :math:`T_\text{max}` with 30 samples
between :math:`\varepsilon_\text{min}` and :math:`\varepsilon_\text{max} = \varepsilon_\text{owc}(T_\text{max})`,
of linear, quadratic, or exponential shape; for :math:`x \in [0, 1]` along the curve the
tension is :math:`T_\text{max}x`, :math:`T_\text{max}x\,(p_2 x + 1 - p_2)`, or
:math:`T_\text{max}(1 - e^{p_2 x})/(1 - e^{p_2})`, with
:math:`\varepsilon_\text{min} = \varepsilon_0 + p_1(\varepsilon_\text{max} - \varepsilon_0)`
(or :math:`\varepsilon_\text{max} - T_\text{max}/p_1` for the linear form with :math:`p_1 \ge 1`).
The shape parameters are the settings-file keys, :math:`p_1 = k_1` and :math:`p_2 = k_2`. The
accepted ranges are :math:`k_1 \ge 0` for the linear shape, :math:`0 \le k_1 < 1` and
:math:`0 < k_2 \le 1` for the quadratic shape, and :math:`0 \le k_1 < 1` and
:math:`10^{-6} \le k_2 \le 700` for the exponential shape (beyond 700 the exponential
overflows; below :math:`10^{-6}` the ratio loses its digits to cancellation).

Each element carries two states: the slow strain :math:`\varepsilon_\text{slow}` and
:math:`T_\text{max}`. With :math:`T_\text{mean} = T_\text{curve}(\varepsilon_\text{slow})`,
:math:`K_1 = \alpha + \beta\,T_\text{mean}`, and ``BA = BA_s|BA_d``,

.. math::

   \dot\varepsilon_\text{slow} =
   \frac{K_1\,\bigl(\varepsilon - \varepsilon_\text{curve}(T_\text{mean})\bigr) + BA_d\,\dot\varepsilon}
        {BA_s + BA_d},
   \qquad
   T = T_\text{mean} + BA_s\,\dot\varepsilon_\text{slow}.

A compressed segment (:math:`\varepsilon < 0`) carries no mean tension, only
:math:`BA_s\,\dot\varepsilon_\text{slow}`, and a tension-only line clamps a negative result
to zero.

The slow strain is eliminated implicitly in every in-step evaluation: a projected, bracketed
backward-Euler solve,
:math:`\varepsilon_\text{slow}^{n+1} = \varepsilon_\text{slow}^{n} + \Delta t\,\dot\varepsilon_\text{slow}(\varepsilon_\text{slow}^{n+1}, \varepsilon, \dot\varepsilon)`
constrained to :math:`\varepsilon_\text{slow} \ge 0`, gives the slow strain that the force
uses. The tangent carries :math:`\partial\varepsilon_\text{slow}^{n+1}/\partial(\varepsilon, \dot\varepsilon)`
from the implicit-function theorem. The step commits that same slow strain, so the force applied
to the nodes equals the reported tension. The working curve is the one regenerated at the
committed :math:`T_\text{max}`. The branch follows the evaluated slow strain: the working curve
applies below the OWC slow strain at :math:`T_\text{max}`, and the OWC applies at or above it,
so the recovered mean tension is continuous across the switch. :math:`T_\text{max}` then
ratchets and regenerates the working curve when the mean tension exceeds it.

**Initial condition.** The static solve uses the rest curve of the initial load history: the
OWC for a virgin rope, or, with ``SYROPE IC``, the working curve regenerated at
:math:`T_\text{max,0}` (zero tension below its zero-tension strain, the OWC above its top
strain :math:`\varepsilon_\text{OWC}(T_\text{max,0})`). A single secant axial stiffness of that
curve is iterated to the solved mean strain. Each element then starts with the running maximum
:math:`\max(T_\text{max,0}, T_e)` (:math:`T_e` for a virgin rope), where :math:`T_e` is its
static tension, and with the slow strain whose kernel tension at the static strain and zero
strain rate equals :math:`T_e`. The dynamic state therefore starts in the static equilibrium,
with no initial transient. Where the static strain lies on the rest curve, that slow strain is
the curve's rest point. Elsewhere, for example along a sagging line, it creeps at
:math:`(T_e - T_\text{mean})/BA_s`. With fixed line ends the mean tension follows from the
geometry and :math:`T_\text{max,0}`, so :math:`T_\text{mean,0}` does not enter the state. A
value that differs from the static mean tension by more than 1 % is reported on the console.
A strain below the working curve's zero-tension strain (a slack rope on that history) stops at
initialisation.

The OWC table must cover the motion. An initial strain outside the table, or a state that
leaves it (a slow strain above the last row's, or a total strain above the last table strain on
the OWC branch), stops with an error rather than extrapolating. The slow strain is well defined
only while the fast spring is stiffer than the working curve: its slope must stay below
:math:`\alpha + \beta T` along the whole curve. For a given table and constants this holds
only for some running maxima. For example, the shipped EXP example constants
(:math:`\alpha = 1.53\times10^8` N, :math:`\beta = 23.12`, :math:`k_1 = 0.2`, :math:`k_2 = 1.5`)
need :math:`T_\text{max} \gtrsim 7.4\times10^5` N. Initialisation scans the table range. A
pretension or ``Tmax0`` outside the admissible range is reported with the offending tension and
the next admissible :math:`T_\text{max}`.

An optional ``SYROPE IC`` row supplies the prior history
:math:`(T_\text{max,0}, T_\text{mean,0})`; see the initial condition above. The model applies
to single-section ``EI = 0`` taut lines; see :doc:`capabilities` for the supported combinations.

Hydrodynamic loads
------------------

Fluid loads follow the Morison equation (:ref:`Morison et al., 1950 <ref-morison1950>`),
decomposed into components normal and tangential to the local unit tangent :math:`\mathbf{t}`.
Reference areas and coefficient names are listed in :doc:`conventions`; :math:`d` is the
hydrodynamic diameter and :math:`A = \pi d^2/4`.

Drag
~~~~

With the relative velocity :math:`\mathbf{u}_r = \mathbf{u}_f - \mathbf{v}` split into
:math:`\mathbf{u}_{r,t} = (\mathbf{u}_r\cdot\mathbf{t})\,\mathbf{t}` and
:math:`\mathbf{u}_{r,n} = \mathbf{u}_r - \mathbf{u}_{r,t}`, the drag per unit length is

.. math::

   \mathbf{f}_D = \tfrac{1}{2}\rho\,d\,C_{d,n}\,\lvert\mathbf{u}_{r,n}\rvert\,\mathbf{u}_{r,n}
                + \tfrac{1}{2}\rho\,\pi d\,C_{d,t}\,\lvert\mathbf{u}_{r,t}\rvert\,\mathbf{u}_{r,t}.

Its derivatives with respect to position (through :math:`\mathbf{t}`) and velocity enter the
implicit tangent exactly.

Added mass and fluid inertia
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The added-mass matrix per unit length and the Froude–Krylov plus fluid-inertia force for a
fluid acceleration :math:`\dot{\mathbf{u}}_f = \dot{\mathbf{u}}_{f,n} + \dot{\mathbf{u}}_{f,t}` are

.. math::

   \mathbf{m}_a = \rho A\left[C_{a,n}\left(\mathbf{I} - \mathbf{t}\mathbf{t}^{\mathsf T}\right)
                + C_{a,t}\,\mathbf{t}\mathbf{t}^{\mathsf T}\right],
   \qquad
   \mathbf{f}_I = \rho A\left[(1 + C_{a,n})\,\dot{\mathbf{u}}_{f,n} + (1 + C_{a,t})\,\dot{\mathbf{u}}_{f,t}\right].

The fluid acceleration is the local (Eulerian) acceleration of the wave field.

Buoyancy and wetting
~~~~~~~~~~~~~~~~~~~~

Self-weight and buoyancy enter as the submerged weight per unstretched metre,
:math:`w = (\rho_A - \rho A)\,g` (positive downward). An ``EQUIVALENT BUOYANCY`` row specifies
:math:`w` and :math:`d` directly and is converted to :math:`\rho_A = w/g + \rho A`.

On the ``EI = 0`` path, an element crossing the free surface uses its wetted fraction
:math:`f`, the linear-interpolation fraction of the chord below the local surface elevation.
Drag, added mass, and fluid inertia act on the wetted part (drag with the wet-weighted relative
velocity), and a recovery load :math:`-(1 - f)\,\tfrac{1}{2}\rho A g L_0` on each end node removes
the buoyancy of the dry part; the derivatives of :math:`f` with respect to the nodal heights are
included in the tangent. On the Hermite path the hydrodynamic loads are integrated over the
deformed element length at three Gauss points, and a Gauss point above the local surface
carries no load. This is exact for fully submerged elements, the regime of lazy-wave cables and
moorings, and approximates the wetted length of a surface-piercing element by the Gauss
sampling. The Hermite buoyancy recovery is exact: the waterline crossings of the cubic
centreline height are located on each element, and the load :math:`\rho g A` is integrated in
closed form against the Hermite shape functions over the dry intervals,

.. math::

   \mathbf{f}_{z} = \rho g A \int_\text{dry} \mathbf{N}^{\mathsf T}\,\mathrm{d}s,
   \qquad
   \frac{\partial \mathbf{f}_{z}}{\partial \mathbf{q}_{z}} = \rho g A\,L_0 \sum_{c}
   \frac{\mathbf{N}(\xi_c)\,\mathbf{N}(\xi_c)^{\mathsf T}}{\lvert z'(\xi_c)\rvert},

where :math:`\xi_c` are the crossings and :math:`z' = \mathrm{d}z/\mathrm{d}\xi`; the tangent is
the exact rate at which the crossings move. The surface is the still waterline of the drag
configuration, or the element's sampled elevation with a host-held field. The Hermite static
solve uses the same recovery at the deck waterline, :math:`z = 0`.

The EI = 0 path assembles the added mass at the current configuration, with its
configuration derivative in the tangent. The Hermite path evaluates the consistent added mass
:math:`\int \mathbf{N}^{\mathsf T}\mathbf{m}_a\mathbf{N}\,\mathrm{d}s` at the start of each step and
holds it over the step. The spatial gradient of the wave field is omitted from the tangent on
both paths; the residual is exact, so the converged state is unaffected.

Waves
~~~~~

Regular waves use linear (Airy) theory (:ref:`Dean & Dalrymple, 1991 <ref-dean1991>`). For
height :math:`H`, period :math:`T`, water depth :math:`h`, and direction :math:`\beta` (the
direction of travel, see :doc:`conventions`), the wavenumber solves the dispersion relation

.. math::

   \omega^2 = g\,k\tanh(k h), \qquad \omega = 2\pi/T,

by monotone bisection in :math:`kh` to a relative interval of :math:`10^{-15}`. With
:math:`\theta = k(x\cos\beta + y\sin\beta) - \omega t`,

.. math::

   \eta = \tfrac{H}{2}\cos\theta, \qquad
   u_h = \tfrac{H}{2}\,\omega\,\frac{\cosh k(z' + h)}{\sinh kh}\cos\theta, \qquad
   u_z = \tfrac{H}{2}\,\omega\,\frac{\sinh k(z' + h)}{\sinh kh}\sin\theta,

where :math:`u_h` acts along :math:`(\cos\beta, \sin\beta, 0)`. The accelerations are the time
derivatives at a fixed point, :math:`\dot u_h = \tfrac{H}{2}\omega^2 \frac{\cosh k(z'+h)}{\sinh kh}\sin\theta`
and :math:`\dot u_z = -\tfrac{H}{2}\omega^2 \frac{\sinh k(z'+h)}{\sinh kh}\cos\theta`. The
vertical velocity is exactly zero at the bed.

**Wheeler stretching** (:ref:`Wheeler, 1970 <ref-wheeler1970>`) evaluates the kinematics at the
stretched coordinate

.. math::

   z' = (z - \eta)\,\frac{h}{h + \eta},

which maps :math:`[-h, \eta]` onto :math:`[-h, 0]`; points above the instantaneous surface
receive zero kinematics. The mapping is undefined when the trough reaches the bed
(:math:`h + \eta \le 0`); that state fails closed with a named error. Without stretching,
:math:`z' = z`.

**Irregular waves** are long-crested sums of Airy components,
:math:`\eta = \sum_i a_i\cos\theta_i`, with :math:`\theta_i = k_i(x\cos\beta + y\sin\beta) - \omega_i t + \phi_i`.
The total elevation is summed first and one Wheeler mapping against that total is applied to
every component. The built-in JONSWAP spectrum (:ref:`Hasselmann et al., 1973 <ref-hasselmann1973>`)
is discretised as follows:

- 200 equal frequency bins on :math:`[0.2\,\omega_p,\ 5\,\omega_p]`, with
  :math:`\omega_p = 2\pi/T_p` and bin width :math:`\Delta\omega = 4.8\,\omega_p/200`; component
  :math:`i` sits at :math:`\omega_i = 0.2\,\omega_p + (i - 1 + u_i)\,\Delta\omega` with
  :math:`u_i` uniform on :math:`(0, 1)`;
- spectral shape
  :math:`S(\omega) \propto \omega^{-5}\exp\!\left[-\tfrac{5}{4}(\omega_p/\omega)^4\right]\gamma^{r}`,
  :math:`r = \exp\!\left[-(\omega/\omega_p - 1)^2/(2\sigma^2)\right]`, with
  :math:`\sigma = 0.07` for :math:`\omega \le \omega_p` and :math:`0.09` above, :math:`\gamma \ge 1`;
- amplitudes :math:`a_i = \sqrt{2 S(\omega_i)\,\Delta\omega}`, with :math:`S` scaled so that the
  discrete zeroth moment is :math:`m_0 = H_s^2/16`;
- phases :math:`\phi_i = 2\pi v_i` with :math:`v_i` uniform on :math:`(0, 1)`.

The draws :math:`u_i, v_i` come from the ``WaveSeed`` option (default 1) through the
Park–Miller MINSTD generator (multiplier 48271, modulus :math:`2^{31} - 1`), so a given
:math:`(H_s, T_p, \gamma, \beta)` and seed always produce the same sea on every platform, and
different seeds give independent realisations. Because the component frequencies are not
commensurate, the record does not repeat (an equally spaced comb would repeat after
:math:`2\pi/\Delta\omega \approx 42\,T_p`). To reproduce a specific external realisation, supply
its component table (amplitude, frequency, phase), for example a WaterKin ``WaveKinMod`` 1
elevation history, which is reduced once to its Fourier components on the run time base.

**Regular nonlinear waves.** A ``stream`` wave is Dean's stream-function wave
(:ref:`Dean, 1965 <ref-dean1965>`) solved numerically by the Fourier method of
:ref:`Rienecker & Fenton (1981) <ref-rienecker1981>`. In the frame moving with the celerity
:math:`c`, with :math:`X = x' - ct` and :math:`Y = z + h` measured from the bed,

.. math::

   \psi(X, Y) = -\bar U Y + \sum_{j=1}^{N} B_j\,\frac{\sinh jkY}{\cosh jkh}\cos jkX,

and :math:`k`, the surface elevations :math:`\eta_m` at :math:`X_m = m\pi/(Nk)`
(:math:`m = 0 \ldots N`, crest to trough), :math:`B_j`, :math:`\bar U`, :math:`Q` and
:math:`R` satisfy the kinematic condition :math:`\psi(X_m, h + \eta_m) = -Q` and the Bernoulli
condition :math:`\tfrac12\lvert\nabla\psi\rvert^2 + g\eta_m = R` on the surface, a zero mean
level, :math:`\eta_0 - \eta_N = H`, and :math:`k\bar U T = 2\pi` (no Eulerian mean current,
:math:`c = \bar U`). Newton's method solves the :math:`2N + 5` equations in units of
:math:`h` and :math:`g`, stepping the height up from a linear wave; ``StreamOrder`` sets
:math:`N` (default 20). The fixed-frame velocity is
:math:`u = \sum jkB_j \cosh(jkY)/\cosh(jkh)\cos jkX`,
:math:`w = \sum jkB_j \sinh(jkY)/\cosh(jkh)\sin jkX`, the local acceleration
:math:`-c\,\partial/\partial X` of these, the surface elevation the cosine interpolant of the
:math:`\eta_m`, and the dynamic pressure :math:`p/\rho + gz = R - \tfrac12\lvert\nabla\psi\rvert^2`.
The field is exact up to the free surface, so no stretching is applied; points above the
surface see still water, and the start-up ramp scales the whole field. Airy and irregular
seas keep the Wheeler stretching above.

**Other spectra, spreading and multiple trains.** Besides JONSWAP, a spectral train can be
ISSC / Pierson–Moskowitz (Bretschneider, :ref:`Pierson & Moskowitz, 1964 <ref-pierson1964>`),
Torsethaugen two-peak (:ref:`Torsethaugen & Haver, 2004 <ref-torsethaugen2004>`), or
Ochi–Hubble two-peak (:ref:`Ochi & Hubble, 1976 <ref-ochi1976>`). With
:math:`\omega_p = 2\pi/T_p`, the one-sided densities are

.. math::

   S_{PM}(\omega) = \tfrac{5}{16} H_s^2 \omega_p^4 \omega^{-5}
   \exp\!\left[-\tfrac{5}{4}(\omega_p/\omega)^4\right],\qquad
   S_{J} = (1 - 0.287\ln\gamma)\,S_{PM}\,\gamma^{r},

   S_{OH}(\omega) = \frac{1}{4}\sum_{j=1}^{2}
   \frac{\left(\tfrac{4\lambda_j + 1}{4}\omega_{pj}^4\right)^{\lambda_j}}{\Gamma(\lambda_j)}
   \frac{H_{sj}^2}{\omega^{4\lambda_j + 1}}
   \exp\!\left[-\tfrac{4\lambda_j + 1}{4}(\omega_{pj}/\omega)^4\right].

The Torsethaugen sea follows DNV-RP-C205 (:ref:`DNV, 2021 <ref-dnvc205>`): with
:math:`T_f = 6.6\,H_s^{1/3}`, a wind-dominated sea (:math:`T_p \le T_f`) has the primary
wind sea :math:`H_1 = r H_s`, :math:`r = 0.7 + 0.3\exp[-(2\epsilon_l)^2]`,
:math:`\epsilon_l = (T_f - T_p)/(T_f - 2\sqrt{H_s})`,
:math:`\gamma_1 = 35\,(2\pi H_1/g T_p^2)^{0.857}`, and a swell of
:math:`H_2 = \sqrt{1 - r^2}\,H_s` at :math:`T_{p2} = T_f + 2` s; a swell-dominated sea has
:math:`r = 0.6 + 0.4\exp[-(\epsilon_u/0.3)^2]`, :math:`\epsilon_u = (T_p - T_f)/(25 - T_f)`,
:math:`\gamma_1 = 35\,(2\pi H_s/g T_f^2)^{0.857}(1 + 6\epsilon_u)` and a wind sea at
:math:`T_{p2} = 6.6\,H_2^{1/3}` (:math:`\epsilon` clamped to [0, 1], :math:`\gamma_1 \ge 1`,
:math:`\gamma_2 = 1`). Each partition is
:math:`E_j G_0 A_{\gamma} \tilde\omega^{-4}\exp(-\tilde\omega^{-4})\gamma_j^{r}` with
:math:`\tilde\omega = \omega T_{pj}/2\pi`, :math:`E_j = H_j^2 T_{pj}/(32\pi)`,
:math:`G_0 = 3.26` and :math:`A_\gamma = (1 + 1.1(\ln\gamma)^{1.19})/\gamma`.

A spectral train is synthesised like the JONSWAP sea: ``WaveComponents`` (default 200) jittered
components over :math:`[0.2\,\omega_{p,\min}, 5\,\omega_{p,\max}]` of its partition peaks,
random phases, and amplitudes scaled so the train's discrete :math:`m_0` is exactly
:math:`H_s^2/16` (:math:`(H_{s1}^2 + H_{s2}^2)/16` for Ochi–Hubble). A train with a spreading
exponent :math:`s > 0` spreads its energy with the cos-2s function

.. math::

   D(\theta) = K(s)\cos^{2s}(\theta - \theta_p),\quad |\theta - \theta_p| \le \pi/2,\qquad
   K(s) = \frac{\Gamma(s + 1)}{\sqrt{\pi}\,\Gamma(s + 1/2)},

discretised into ``WaveDirections`` (default 9) equal-angle bins whose weights are the bin
integrals of :math:`D`; each direction has its own jittered frequency set, and component
:math:`i` travels along its own heading :math:`\beta_i`. Several trains (``wavetrain`` rows) are
superposed component by component, each with its own heading and seed offset, and the
kinematics of the whole sea follow the irregular-wave conventions above: the total elevation
first, one Wheeler mapping against it, then the component sums. A long-crested JONSWAP train
reproduces the JONSWAP synthesis above for the same seed.

The optional ``rampTime`` :math:`T_r` multiplies every wave amplitude by
:math:`r(t) = \tfrac{1}{2}\left[1 - \cos(\pi t/T_r)\right]` for :math:`t < T_r`; the fluid
acceleration then includes the ramp rate, :math:`r\,\dot{\mathbf u} + \dot r\,\mathbf u`. The
current is not ramped, because the static initial condition already includes it.

Currents
~~~~~~~~

A current is uniform or given as a profile of :math:`N \ge 2` levels
:math:`(z_j, \mathbf{u}_j)` with strictly increasing :math:`z_j`. The velocity is interpolated
linearly between levels and held at the end values outside the profile. A WaterKin
``CurrentMod`` 1 depth table uses the same interpolation. In OpenFAST (maintained by NLR, the National Laboratory of the
Rockies, formerly NREL), the SeaState field supplies waves and current; SeaState's standard ``CurrMod = 1`` steady current (sub-surface
power law :math:`U_{ss}\,[(z + h)/h]^{1/7}`, linear near-surface profile, and depth-independent
part) is reproduced so that wave and current contributions can be separated when a WaterKin
file selects them independently. ``CurrMod = 2`` user currents cannot be separated and are
rejected whenever separation is required.

The fluid velocity seen by the drag is the sum of the current and the wave velocity. In a
coupled run the fluid velocity, acceleration, and local surface elevation are sampled by the
host at the nodes and held over the step.

Seabed contact
--------------

Contact acts at the line nodes on the translational DOFs. With the penetration
:math:`g = z_\text{floor}(x, y) - z` and the nodal stiffness :math:`k_n` (:doc:`conventions`),
the normal force is a :math:`C^1` penalty law with blend width :math:`\delta = 10^{-6}\,\mathrm{m}`:

.. math::

   f_n^{(k)}(g) =
   \begin{cases}
     0, & g \le -\delta, \\[2pt]
     k_n\,\dfrac{(g + \delta)^2}{4\delta}, & -\delta < g < \delta, \\[6pt]
     k_n\,g, & g \ge \delta .
   \end{cases}

The force and its tangent are continuous through first contact, which avoids an active-set
kink in Newton's method; beyond the blend the law is exactly linear. A normal damper acts only
while the node moves into the bed,

.. math::

   f_n^{(c)} = -c_n\,\chi(g)\,\min(v_z, 0),

where :math:`\chi` rises linearly from 0 to 1 across the blend. With the total normal reaction
:math:`N = f_n^{(k)} + f_n^{(c)}`, seabed friction is a stick-slip spring on the horizontal
motion of every node, on the ``EI = 0`` and the finite-EI lines alike. Node :math:`i` holds a
horizontal spring of stiffness :math:`k_i`, the nodal normal contact stiffness (the shear
stiffness equals the normal stiffness, OrcaFlex's default), to its anchor
:math:`\mathbf{a}_i`. With :math:`\mathbf{d} = \mathbf{x}_i - \mathbf{a}_i` the horizontal
stretch, the spring sticks and slides:

.. math::

   \mathbf{f}_t = \begin{cases} -k_i\,\mathbf{d}, & k_i\lvert\mathbf{d}\rvert \le \mu N,\\
   -\mu N\,\mathbf{d}/\lvert\mathbf{d}\rvert, & \text{otherwise},\end{cases}

and at every committed state a sliding spring's anchor moves to the capacity distance behind
the node (return mapping), so the force carries over exactly into the next step; a node that
leaves the seabed takes its anchor with it. The anchors are state: they roll back with a
failed or subdivided step, join the aggregate snapshot, and are part of the OpenFAST
checkpoint. The force and its position and velocity Jacobians (through :math:`N`) are exact.
Unlike a velocity-regularised Coulomb law, the spring holds its force at rest.

In statics without a current the anchors are the static node positions and carry no force. In
a deck current the static solve carries friction springs from a reference position, capped
smoothly at :math:`\mu N`,

.. math::

   \mathbf{f}_t = -\frac{k_i\,\mathbf{d}}{\sqrt{1 + (k_i\lvert\mathbf{d}\rvert/\mu N)^2}},

with :math:`\mathbf{d}` the horizontal displacement from the reference. The reference is the
still-water laid shape: the same line is solved without the current first, as OrcaFlex
references its friction to the laid position, so the equilibrium depends on that path. The
springs are memoryless, so the equilibrium does not depend on how the current is brought in;
on the finite-EI path their capacities are frozen with the drag in the fixed-point passes and
live in the final Newton solve, on the ``EI = 0`` path they are frozen at the line-search base
point of the energy step. The anchors then start at :math:`\mathbf{a}_i = \mathbf{x}_i -
\mathbf{f}_s/k_i` for the static spring force :math:`\mathbf{f}_s`, so the march starts in
the static equilibrium and a line held at rest in the current keeps its static end tensions.
Standalone and coupled runs share this model.

**Anisotropic friction.** With ``frictionMuAxial`` :math:`\mu_a` and ``frictionMuLateral``
:math:`\mu_n` (the OrcaFlex axial and normal coefficients) the capacity depends on the slip
direction. With :math:`\hat{\mathbf a}` the unit horizontal projection of the local line axis
(the nodal tangent on the finite-EI path, the chord through the neighbouring nodes on the
``EI = 0`` path), :math:`\hat{\mathbf n} = \mathbf e_z\times\hat{\mathbf a}`, and
:math:`d_a = \mathbf d\cdot\hat{\mathbf a}`, :math:`d_n = \mathbf d\cdot\hat{\mathbf n}`,

.. math::

   \mu(\mathbf d) = \frac{\sqrt{(\mu_a d_a)^2 + (\mu_n d_n)^2}}{\lvert\mathbf d\rvert},

so a slip along the line is resisted by :math:`\mu_a N`, one across it by :math:`\mu_n N`, and
the force stays collinear with :math:`\mathbf d`. Both laws above use the capacity
:math:`\mu(\mathbf d) N` in place of :math:`\mu N`, the return mapping included. The line axis
is evaluated at every residual evaluation; the Jacobians carry the dependence of the capacity
on the slip direction but omit its dependence on the line axis.
A vertical axis (no horizontal projection) takes
:math:`\mu_n`. An equal pair runs the isotropic law unchanged.

A structured bathymetry grid gives :math:`z_\text{floor}(x, y)` by bilinear interpolation,
clamped to the nearest edge outside the grid. On a sloped patch the elastic contact is
frictionless and normal to the surface, on the ``EI = 0`` and the Hermite path alike: with
:math:`s = \sqrt{1 + z_{,x}^2 + z_{,y}^2}` the force is :math:`f_n^{(k)}(g/s)\,\mathbf{n}`,
:math:`\mathbf{n} = (-z_{,x}, -z_{,y}, 1)/s`, and its tangent
:math:`f_n^{(k)\prime}\,\mathbf{n}\mathbf{n}^{\mathsf T}` (the grid curvature is not
differentiated); on a level floor this is exactly the vertical law above. Its damping acts on
the velocity along :math:`\mathbf{n}`, and a friction capacity is :math:`\mu` times its
magnitude. The static solve and the dynamics use the same law, so a grounded run on a slope
rests at a static equilibrium of the march; its tension changes along the slope by the weight
component :math:`w\sin\theta` per unit length, and a run pushed down the slope into its anchor
is compressed there. On the Hermite path the contact acts on the nodal position DOFs only, not
on the tangent DOFs.

Static equilibrium
------------------

Every run starts from a converged Newton equilibrium of
:math:`\mathbf{R}(\mathbf{q}) = \mathbf{f}_\text{int}(\mathbf{q}) - \mathbf{f}_\text{ext}(\mathbf{q}) = \mathbf{0}`
on the free DOFs, not from dynamic relaxation. External loads are the submerged weight, the
seabed reaction, end-connection
springs, on the ``EI = 0`` path the steady current, and, on the Hermite path, prescribed nodal
loads. The seabed reaction enters every residual and tangent evaluation; only the load
continuation of the ``EI = 0`` path (step 4 below) scales it during its intermediate stages.
Waves act from the start of the dynamic march, whose initial acceleration is computed
consistently with them.

.. _static-ei0:

EI = 0 path
~~~~~~~~~~~

The static-only deck route and the dynamic initial condition (``TMax = 0`` included) call the
same initialiser, so both report the same equilibrium. An ``EI = 0`` line cannot carry
compression: the accepted state is the tension-only equilibrium, with no element compressed
beyond a :math:`10^{-6}` relative round-off band; a slack element reports zero tension.

1. **Seeds**, tried in order until one gives an admissible equilibrium: the analytical
   extensible catenary (grounded when the anchor lies within
   :math:`\max(1\,\mathrm{mm}, 10^{-5}\lvert z_\text{floor}\rvert)` of the seabed, on the seabed
   plane through the anchor with its slope along the span; otherwise fully suspended); on a
   bathymetry, a grounded run draped on the actual seabed profile that lifts off tangentially
   into a suspended catenary; the catenary on the seabed lowered by the static penalty
   penetration :math:`w/(k_\text{Bot} d)` (a soft seabed); for a grounded line longer than span
   plus rise, the :math:`H \to 0` limit (a straight hanging leg with the excess on the seabed:
   where the seabed keeps falling past the fairlead the excess slides down the frictionless
   slope and folds back, with zero tension at the fold; otherwise it lies slack); the other
   catenary variant; and a two-leg polyline of the line length that never compresses the line.
   Seed nodes deeper below the seabed than their static penetration are lifted onto it.
2. **Compression-capable Newton** at the full load (Armijo line search on
   :math:`\tfrac{1}{2}\lVert\mathbf{R}\rVert^2`, :ref:`Armijo, 1966 <ref-armijo1966>`). A converged
   state with every element in tension is the tension-only equilibrium and is accepted; a
   compressed strut or arch is rejected.
3. **Tension-only Newton** otherwise: the compressed branch gets a stiffness :math:`r\,EA`, with
   :math:`r\,EA` of order :math:`10^{-3}` and then :math:`10^{-7}` of the largest nodal load (too
   weak to hold an arch), and finally :math:`r = 0`. Each stage is a Levenberg-Marquardt Newton
   iteration whose line search descends the potential energy (strain energy, seabed penalty
   potential, load potential; a steady-current drag is frozen at the base point of each search),
   so it approaches stable equilibria only; its shift never drops below :math:`10^{-12}` of the
   geometric stiffness scale, which keeps directions without stiffness (a grounded node between
   slack elements) solvable.
4. **Load continuation** in the stages 0.25, 0.5, 0.75, 1.0 is the fallback of each law. Each
   stage scales the seabed stiffness with its load factor, which keeps a stiff contact well
   conditioned from a cold seed; the final stage restores the full stiffness, so the converged
   equilibrium is unchanged.
5. **Convergence.** :math:`\lVert\mathbf{R}\rVert_\infty / s < 10^{-8}` or
   :math:`\lVert\mathbf{R}\rVert_\infty < 10^{-14}`, where :math:`s` is the largest of the
   infinity norms of the external load, the seabed reaction, and any current drag; a residual
   that stops decreasing inside the round-off floor
   :math:`16\,\epsilon\,\max(EA/L_0)\max(\lvert q\rvert, L_0)` of the nodal forces is converged.
   A line with no admissible equilibrium fails closed, naming the line, the seeds and stages
   tried, and the reason.

A steady current (uniform, or a depth profile sampled at each node's own elevation as in the
march) enters the static equilibrium as Morison drag on the line at rest, so the out-of-plane
offset is part of the initial condition. ``Connect`` and ``Free`` points are then placed at
their static force balance: a trust-region Newton iteration on the point positions whose
residual sums the attached lines' end forces (each line re-solved for the moving ends), the
point's weight and buoyancy, and its current drag. It is the rigid-body equilibrium solve
described under *Rigid-body initial condition* below, with points as the only unknowns.

Hermite path
~~~~~~~~~~~~

A finite-bending deck cable is initialised by the practice of "catenary statics, then full
statics" (``continuation cable_statics``, the default; ``sequenced cable_statics`` runs the
mesh-sequenced route described below first and this route as its fallback):

1. **Exact EI = 0 seed.** The multi-segment extensible catenary of the same sections (signed
   submerged weight, so buoyancy sections are included) is solved through the endpoints, on the
   declared seabed plane (flat ``WtrDpth``, or the least-squares plane of a structured bathymetry
   along the chord) or hung in mid-water when there is no seabed. The solve runs in the vertical
   plane of the chord, translated to the anchor, so a translated or rotated deck gives the same
   numbers to round-off. When the only catenary through the endpoints passes below a planar
   seabed (flat, or a structured bathymetry whose profile under the chord is straight), the
   line is too long for its span to hang as a catenary. Without buoyant sections the line is a
   convex curve between the chord and the seabed, so it cannot be longer than the path down to
   the seabed, along it to below the fairlead and up; a longer line is in the same case. Its
   bending stiffness can still hold an equilibrium: the excess length stands up as a bow in
   the vertical plane, carried in compression by the bending stiffness. That state is reached
   by displacement control, as a line is laid: the anchor is moved back along the chord (by the
   excess length plus a tenth of the line, then half as far again until the catenary exists),
   the line is solved there by the steps below, and the anchor is walked back to its position
   in steps, each a minimisation of the potential energy from a secant predictor (plain Newton
   in a current, whose drag is not a potential force), the step halved on a failure and
   doubled after a success. The walked state is polished, audited and checked for stability
   like any other; the driver then reports its element-mean compression. When the walk ends
   short of the anchor (the bow snaps through, or the line would have to lie on itself), or the
   bathymetry does not continue its plane behind the anchor, the run stops and names the
   geometry: no admissible planar equilibrium clear of the seabed was found, the line being too
   long for its span to hang as a catenary, with where the walk stopped. A walked state that
   the audit below rejects (self-contact, instability) hands over to the mesh-sequenced route,
   and is named the same way when that finds none; with ``tensile_safety True`` a walked state
   in compression is rejected by name (a bow held by bending is compressed).
   The walked arch faces the anchor. When the walk fails or its state is rejected, a line in
   still water over a flat seabed is also traced to the bow on the other side of the line
   hanging from the fairlead, which the walk cannot reach without passing the arch through
   that line. The anchor is placed beyond the fairlead, as far from it as the walk moved it
   back, where the arch faces away from the real anchor, and walked toward the fairlead. Each
   step is bounded so as not to pass the point where the grounded run, folded back about its
   touchdown point, reaches the real anchor (to 2 % of the span), and the walk stops there or
   at the last state before the arch turns over. The grounded run is folded, the corner at
   the touchdown relaxed by an energy minimisation, and the anchor walked the rest of the way
   into place; the state is polished, audited and checked for stability like the walked one.
   The bow is compressed, so it is not traced with ``tensile_safety True``, nor in a current.
   On a mesh with elements shorter than :math:`\sqrt{EI/EA}` the continuation below has an
   iteration budget of :math:`2\times10^{6}` divided by the element count (a solve that
   converges there does so within it), and no other route follows a failure: the run stops in
   seconds and names the mesh.
2. **Continuation in** ``EI``. The family :math:`\lambda\,EI` is traced from a fraction whose
   bending length :math:`(\lambda EI/|w|)^{1/3}` is below a quarter of the element length up to
   :math:`\lambda = 1`, every load at full value. Each stage minimises the total potential energy
   (elastic, weight, buoyancy lost above the surface, seabed penalty, end springs) by a modified
   Newton iteration: the tangent is shifted until a banded Cholesky factorisation succeeds, so
   every direction descends toward a *stable* equilibrium, and the step length follows an Armijo
   line search on the energy. Close to equilibrium full Newton steps on the complete tangent
   take over, because the seabed contact and waterline laws are only piecewise smooth (and the
   curvature of a structured bathymetry grid is not differentiated). Stages that jump
   further than a small fraction of the line are retried at a smaller step; a first stage that
   does not converge is restarted at full ``EI``.
3. **Polish and diagnosis** on the deck mesh at the static tolerance. A residual floor above it
   is accepted at :math:`10^{-6}` (or, for a line on the seabed, :math:`10^{-5}`, nodes
   alternating across the 1 µm contact blend) and reported. A deck mesh too coarse for the
   bending length (continuation or polish failing, or an element folded beyond
   :math:`L_0\kappa = 0.7`) repeats the route on a mesh refined until the elements are at most
   half the bending length (when the continuation on that mesh runs out of its iteration
   budget, on the coarser powers of two down to twice the deck mesh); a rotationally
   restrained end whose boundary layer :math:`\sqrt{EI/T}` is shorter than its element, which
   shows as a polish that does not converge or a compression at that end, repeats it on a mesh
   graded geometrically toward that end. The refined mesh is reported.
4. **Physical-branch audit.** The converged state must be clear of self-contact (no two
   non-adjacent elements closer than their diameters), unfolded, and stable (positive definite
   tangent). Element-mean axial compression of a stable bent equilibrium is admissible and is
   reported; with ``tensile_safety True`` it is rejected. That contract is judged, as in the
   dynamics, on the element-mean axial force averaged over each element and its two
   neighbours, compared with the tighter of the strain band :math:`2\times10^{-6}\,EA` (with the
   ``EA`` of the centre element) and a thousandth of the peak tension (the largest pointwise
   resultant); the pointwise resultant of a coarse element oscillates about that mean (tens of
   kN on the Gulf of Mexico 80 m lazy wave at 48 elements) without entering equilibrium, and
   does not drive the mesh refinement. Before rejecting, the static solve
   splits the elements whose three-element mean falls below the band (and their neighbours) and
   solves again, since a coarse element can dip where the line is in tension; it stops, and the
   solve names the compression, when two passes do not relieve it, when it is the reaction of a
   held end (the minimum settled on the end node over two passes, which no splitting changes),
   when the polish on the split mesh does not converge, or at an element budget. The report is
   judged on element-mean forces: the mean over each element and its two neighbours
   (length-weighted) beyond the
   admissible band, or a node tension (the segment tension, :doc:`outputs`) in compression by
   more than a tenth of the peak tension.
   The nodal seabed-contact reaction kinks the line at the touchdown node, and the
   element there can carry a small negative mean between tensile neighbours at any mesh; that dip
   does not enter equilibrium and is not reported. A failing route falls back to the
   mesh-sequenced solve below, audited alike. On the Lozon cables the two routes give the same
   equilibrium to the solver tolerance (positions within 2e-5 m, curvature within 1e-6 1/m).

A deck current loads the line at rest with its steady Morison drag (the relative velocity is
the current itself; the element fluid velocity is the mean of the current at its two nodes,
as in the dynamics). The drag follows the line and is not a potential force: it is carried
along the stiffness continuation frozen at the last accepted configuration, then converged
at full stiffness by a fixed-point iteration on the frozen drag and a final Newton solve with
the configuration-dependent drag and its Jacobian. If the frozen-drag continuation fails, the
still-water equilibrium is found first and the drag is brought in by load steps. A load step
that stalls hands over to a pseudo-arclength continuation in the drag fraction
:math:`\lambda`: each corrector solves the bordered system
:math:`[K,\,-D;\ t_q^T,\ t_\lambda]` (:math:`D` the drag at full current) by block elimination
with two banded solves, so the path can turn at a limit point in :math:`\lambda`. A current
across the vertical plane of the line makes the static solve three-dimensional, as the
dynamics are.

Without seabed friction, a current that pushes the grounded run toward its anchored end
harder than the suspended span pulls it back puts that run into compression; a cable cannot
hold it, and the run would have to buckle or slide past its end. There is no tension-only
equilibrium that continues the still-water layout in its plane. With the current in the
vertical plane of the line the problem is symmetric and the planar branch ends in a
pitchfork: the equilibrium beyond it turns the grounded run sideways on the seabed, past its
anchor and back (a lazy wave in a 0.7 to 1.2 m/s current along its plane, for example).
Newton cannot leave the plane of a symmetric problem, so the line is then solved in 3D under
the current with a cross-flow of 1 % of its horizontal speed, and that state is polished at
the exact current and audited like any other. A planar equilibrium in such a current is also
checked with the line free to leave its plane; one that is a saddle there (the grounded run
held in the plane only by the symmetry) is replaced by this 3D solve; when the 3D solve finds
no equilibrium, the planar one is kept and a note says that it is held in the plane of the line
and unstable out of it. Where the planar path
ends depends on the element size (whether a load step or a polish crosses the limit), so the
branch is selected on the deck mesh coarsened by halving every section's ``NumSegs`` as long as
the elements stay within one bending length :math:`(EI/|w|)^{1/3}`: that coarse line goes
through every route above, and its equilibrium is carried to the deck mesh by cubic-Hermite
interpolation along the unstretched arc length and polished there at the exact current, free
to leave the plane, then audited and checked for stability. Decks whose meshes are halvings of
one another select the same coarse mesh, and so the same equilibrium; a residual floor above
the tight tolerance on a fine mesh is accepted at the coupled static tolerance and noted. The
result can be a folded (hairpin) layout:
a very strong current on a frictionless seabed carries the line past its anchor and it folds
back to it, either as a loop lying on the seabed with its two legs side by side or, for a
lazy wave, with the buoyant arch carried beyond the anchor and the line coming back down to it.
Such a layout is a genuine equilibrium: the drag on the legs is held by the tension of both,
the seabed carries the grounded part, and the tangent is positive definite. OrcaFlex statics
seeded with the CableDyn shape converge to the same layout, and a dynamic run from it stays at
rest and returns to it after a disturbance. The standalone driver prints a note when the
grounded run of the static layout lies past its anchor (``the static layout carries the
grounded run ... m past its anchor``), because the layout differs from the usual one and a
real seabed with friction, or a shorter line, would not reach it. When the 3D solve also
fails, the solve stops with
the diagnosis (the element-mean compression, where it sits, and its height above the seabed)
and suggests ``frictionMu``; neither the pseudo-arclength continuation of the planar path
nor the mesh-sequenced route, which reach the same planar limit, is then run. With friction the grounded run is held; a lazy wave whose buoyant arch the
current pushes toward the touchdown then goes into compression at the touchdown instead,
which friction cannot hold, and the solve names that. A line that is already compressed
without the current is named as such. A lazy wave whose seabed run points downstream reaches
these limits at a fraction of the current that the same layout holds when the current comes
from the other side.

The free surface uses the partial immersion of the circular section: over the band
:math:`|z| < D/2` the displaced area, and so the buoyancy lost above the water, follows the
circular-segment fraction of the section, which gives a line floating at the surface its
waterplane stiffness. The static and dynamic paths share this law.

The fallback and the dynamic Hermite path use

- **continuation in** ``EI``, and, for a net-buoyant line, **continuation in buoyancy**, which
  tracks the smooth branch of the buoyant arch from the seed;
- the **consistent Hermite self-weight**, :math:`wL_0/2` on each end :math:`r_z` and
  :math:`\pm wL_0^2/12` on the end :math:`m_z`;
- **mesh sequencing**: the equilibrium is selected on a coarse mesh that preserves section,
  load, and contact-diameter interfaces, prolonged through the Hermite interpolant, and polished
  on each finer level and finally on the exact user mesh with the same residual test. Deterministic
  fallback hierarchies are tried before a direct cold solve; every attempt uses the same
  physical inputs and tolerance;
- a **dimensionless per-DOF residual norm**, with tangent-DOF residuals divided by the nodal
  tributary length.

The static tolerance is :math:`\min(10^{-6},\ \text{tol}_\text{dyn})`, where
:math:`\text{tol}_\text{dyn}` is the dynamic Newton tolerance, so the dynamics start from an
equilibrium at least as tight as their own stopping criterion. Two checks keep a numerically
balanced but unresolved state out of the dynamics: an element whose sampled
:math:`\max(L_0\,\kappa) > 0.7` is an unresolved fold and is rejected, and, with
``adaptive_mesh``, exhausting the refinement cap while the curvature or contact topology remains
under-resolved is reported.

Coupled systems of lines, points, and bodies are initialised in sequence: non-line objects are
held while each line is solved, then released and the coupled system is solved.

Modal analysis
--------------

With ``nModes`` each line's natural frequencies and mode shapes about its static equilibrium
solve the generalised symmetric eigenproblem

.. math::

   \mathbf K\,\boldsymbol\phi = \omega^2\,\mathbf M\,\boldsymbol\phi

over the free degrees of freedom, the two ends held (a fixed anchor and a held or prescribed
fairlead). :math:`\mathbf K` is the static tangent stiffness: for an ``EI = 0`` element of rest
length :math:`l_0`, current length :math:`l` along :math:`\mathbf e`, and tension
:math:`T = EA(l/l_0 - 1)` (zero on a slack tension-only element),
:math:`\mathbf k = (EA/l_0)\,\mathbf e\mathbf e^{\mathsf T} + (T/l)(\mathbf I - \mathbf e\mathbf e^{\mathsf T})`
assembled as :math:`[[\mathbf k, -\mathbf k], [-\mathbf k, \mathbf k]]`; for a cubic-Hermite
element, the Hessian of its strain energy. The linearised seabed normal contact adds its
stiffness on the contacting nodes. :math:`\mathbf M` is the consistent structural mass plus the
Morison added mass at the static shape. Friction, damping, and drag are left out, as in any
undamped modal analysis. :math:`\mathbf K` and :math:`\mathbf M` are banded (the nodes are
numbered along the line), and only the ``nModes`` lowest modes are computed: LAPACK ``dsbgvx``
gives the eigenvalues, inverse iteration with the banded factors of
:math:`\mathbf K - \omega^2\mathbf M` gives each shape (M-orthogonal to the shapes before
it, so a repeated frequency gets independent shapes), and the returned :math:`\omega^2` is the
Rayleigh quotient of the shape, summed in quadruple precision. That quotient keeps the lowest
frequencies of a long, axially stiff line accurate, where the band reduction alone loses
digits. Memory grows linearly with the line length. The frequencies
:math:`f = \omega/2\pi` and the nodal translations of each mode (scaled to a unit largest
nodal displacement) are written to ``<root>.modes.out``. For a taut string of tension
:math:`T` the transverse frequencies are :math:`f_n = (n/2L)\sqrt{T/m}`; for a shallow sagging
cable the out-of-plane and antisymmetric in-plane modes follow the string with the horizontal
tension :math:`H`, and the symmetric in-plane modes satisfy
:math:`\tan(\beta/2) = \beta/2 - (4/\lambda^2)(\beta/2)^3`, :math:`\omega = (\beta/L)\sqrt{H/m}`,
with :math:`\lambda^2 = (wL/H)^2\,L/(H L_e/EA)` (:ref:`Irvine & Caughey, 1974 <ref-irvine1974>`);
the modal analysis reproduces both.

Time integration
----------------

Both paths advance

.. math::

   \mathbf{M}\ddot{\mathbf{q}} + \mathbf{f}_\text{int}(\mathbf{q}, \dot{\mathbf{q}})
   = \mathbf{f}_\text{ext}(\mathbf{q}, \dot{\mathbf{q}}, t)

with the generalised-α method (:ref:`Chung & Hulbert, 1993 <ref-chung1993>`). There are no
rotation variables, so :math:`\mathbf{q}` and :math:`\ddot{\mathbf{q}}` share one vector space and
no gyroscopic or rotation-chart terms arise. The parameters follow from the spectral radius
:math:`\rho_\infty` (``rhoInf``, :doc:`conventions`):

.. math::

   \alpha_m = \frac{2\rho_\infty - 1}{\rho_\infty + 1}, \qquad
   \alpha_f = \frac{\rho_\infty}{\rho_\infty + 1}, \qquad
   \beta = \tfrac{1}{4}(1 - \alpha_m + \alpha_f)^2, \qquad
   \gamma = \tfrac{1}{2} - \alpha_m + \alpha_f .

With :math:`\mathbf{q}_{n+1}` as the Newton unknown,

.. math::

   \mathbf{a}_{n+1} = \frac{\mathbf{q}_{n+1} - \mathbf{q}_n}{\beta\,\Delta t^2}
                    - \frac{\mathbf{v}_n}{\beta\,\Delta t}
                    - \left(\frac{1}{2\beta} - 1\right)\mathbf{a}_n, \qquad
   \mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\left[(1 - \gamma)\,\mathbf{a}_n + \gamma\,\mathbf{a}_{n+1}\right],

and the residual is evaluated at the intermediate states

.. math::

   \mathbf{G} = \mathbf{M}\left[(1 - \alpha_m)\,\mathbf{a}_{n+1} + \alpha_m\,\mathbf{a}_n\right]
              + \mathbf{f}_\text{int}(\mathbf{q}_\alpha, \mathbf{v}_\alpha)
              - \mathbf{f}_\text{ext}(\mathbf{q}_\alpha, \mathbf{v}_\alpha, t_\alpha),
   \qquad (\cdot)_\alpha = (1 - \alpha_f)(\cdot)_{n+1} + \alpha_f(\cdot)_n .

Time-dependent loads (waves) are evaluated at :math:`t_\alpha`. The Hermite path blends the
*forces* instead (``alpha_force_blend``, default ``True``; ``False`` selects the form above),

.. math::

   \mathbf{G} = \mathbf{M}\left[(1 - \alpha_m)\,\mathbf{a}_{n+1} + \alpha_m\,\mathbf{a}_n\right]
              + (1 - \alpha_f)\,\mathbf{r}(\mathbf{q}_{n+1}, \mathbf{v}_{n+1}, t_{n+1})
              + \alpha_f\,\mathbf{r}(\mathbf{q}_n, \mathbf{v}_n, t_n),
   \qquad \mathbf{r} = \mathbf{f}_\text{int} - \mathbf{f}_\text{ext},

with the same accuracy and tangent. Its position and tangent-vector DOFs rotate with the line,
and a blended configuration :math:`\mathbf{q}_\alpha` would shorten the tangents of a rotating
element, a spurious stretch of order :math:`\alpha_f(1 - \alpha_f)(\omega\Delta t)^2 EA` that
raises the mean tension at large time steps. :math:`\mathbf{r}_n` is the previous step's
converged force, reused unless the state or the load configuration changed in between. Since it
enters every later step, the accepted state is taken one Newton correction beyond the tolerance
(the last factorisation reused), so a loose Newton tolerance does not bias the response. Both
forms converge to the same response; the force blend removes the rotation-driven bias (see
``alpha_force_blend`` in :doc:`options`). The effective tangent is

.. math::

   \frac{\partial\mathbf{G}}{\partial\mathbf{q}_{n+1}}
   = \frac{1 - \alpha_m}{\beta\,\Delta t^2}\,\mathbf{M}
   + (1 - \alpha_f)\,\mathbf{K}_q
   + (1 - \alpha_f)\,\frac{\gamma}{\beta\,\Delta t}\,\mathbf{K}_v ,

where :math:`\mathbf{K}_q = \partial(\mathbf{f}_\text{int} - \mathbf{f}_\text{ext})/\partial\mathbf{q}`
and :math:`\mathbf{K}_v = \partial(\mathbf{f}_\text{int} - \mathbf{f}_\text{ext})/\partial\mathbf{v}`
collect the structural, damping, drag, contact, and friction terms, and :math:`\mathbf{M}`
includes the added mass. On the ``EI = 0`` path the added-mass contribution also carries its
configuration derivative.

**Newton iteration.** An ``EI = 0`` step starts from the Newmark predictor
:math:`\mathbf{q}^{(0)}_{n+1} = \mathbf{q}_n + \Delta t\,\mathbf{v}_n + \Delta t^2(\tfrac{1}{2} - \beta)\,\mathbf{a}_n`
(the acceleration guess :math:`\mathbf{a}^{(0)}_{n+1} = \mathbf{0}`). A Hermite step with the
force blend (the default) starts from the constant-acceleration predictor
:math:`\mathbf{q}^{(0)}_{n+1} = \mathbf{q}_n + \Delta t\,\mathbf{v}_n + \tfrac{1}{2}\Delta t^2\,\mathbf{a}_n`
when the committed acceleration balances the committed force :math:`\mathbf{r}_n`,
:math:`\lVert\mathbf{M}\mathbf{a}_n + \mathbf{r}_n\rVert \le \tfrac{1}{2}\lVert\mathbf{M}\mathbf{a}_n\rVert`
over the translational free DOFs, and from the Newmark predictor otherwise (and always with
``False alpha_force_blend``): while an
under-resolved mode rings, a loosely converged step would keep the extrapolated acceleration,
whereas the zero guess damps it. A failed constant-acceleration attempt is re-solved from the
Newmark predictor. The test uses only the committed state, so a restart is exact.
The ``EI = 0`` path uses an Armijo line search on :math:`\tfrac{1}{2}\lVert\mathbf{G}\rVert^2`
and converges when :math:`\lVert\mathbf{G}\rVert_\infty` relative to the larger of the external
load and the initial residual falls below the relative tolerance, or below the absolute
tolerance. Two near-tolerance cases are also accepted. When the iteration budget is spent or
the line search stalls, a residual within twice the relative tolerance (never looser than
:math:`\sqrt{\epsilon}`) is converged. A residual that stops contracting is converged when it is
below :math:`10^{-4}` of the force scale and within the round-off floor
:math:`\epsilon\,\max\lvert q\rvert\,\lVert\mathbf{A}\rVert_\infty` of the effective operator
:math:`\mathbf{A}`, which no Newton step can reduce (a finely meshed stiff line at small time
steps sits there). The Hermite path uses a backtracking line search that halves the step until the
residual decreases; its convergence norm divides each translational residual by a force scale
frozen at the predictor and each tangent-DOF residual additionally by the nodal tributary
length. The force scale is the largest inertial or residual entry, bounded below by 1 N, by the
peak element weight, and by the round-off floor
:math:`16\,\epsilon\,\max(EA/L_0)\,x_\text{ref}/\text{tol}` of the axial forces
(:math:`x_\text{ref}` the largest coordinate magnitude or element length). The
deck defaults are relative tolerance :math:`10^{-8}`, absolute tolerance :math:`10^{-14}`, 30
iterations, and 12 backtracks (option ``dynamic_solver``). Option ``modified_newton`` reuses the
factorised tangent within a step while each iterate reduces the residual to at most a quarter
of the previous one, and rebuilds it otherwise; it does not relax the tolerance.

A standalone Hermite step on the constant-acceleration branch (so only with the force blend)
also starts from the previous step's factorised tangent, under the same quarter-contraction
rule for at most three iterations, and must reach :math:`0.3` of the tolerance while on it (a
reused factor converges linearly, so this keeps the result as close to the converged solution
as a fresh tangent's last quadratic step). If a step on the reused factor still needs as many
factorisations as the last fresh step, the reuse saved nothing: the following steps then
factorise fresh, for 1 step after the first such step, 2 after the second in a row, and so on,
doubling up to 32 steps. A step that fails on the reused factor is re-solved on a fresh one.
Cables run through the coupling aggregate
(OpenFAST, FAST.Farm, and mixed mooring-and-cable decks) factorise every step, because the
reused factor is not part of the checkpoint state.

**Prescribed motion.** Prescribed DOFs (a moving fairlead or hang-off) are held at their
target position during the Newton solve; their prescribed velocity and acceleration enter the
free-DOF residual through the mass and stiffness coupling. At initialisation, and after a
non-stepping ``motionFile`` boundary commit, the initial acceleration solves
:math:`\mathbf{M}_{ff}\,\mathbf{a}_f = \mathbf{f}_f - \mathbf{M}_{fp}\,\mathbf{a}_p`, so an
accelerated fairlead loads the adjacent free DOFs through the off-diagonal consistent mass. The
core does not ramp prescribed motion; ramp it from rest over at least one wave period to avoid a
start-up transient.

**Adaptive step subdivision.** If a step fails to converge, or (Hermite path) converges but
changes the chord of any element (the relative position of its two nodes) by more than 20 % of
its length, turns any nodal tangent by more than 15°, or shrinks a nodal tangent magnitude
:math:`\lvert\mathbf{m}\rvert` below 0.01 (a singular parametrisation), the step-start state is restored and the
same interval is re-advanced as
:math:`n` substeps of :math:`\Delta t/n`, with :math:`n = 4, 16, 64, 256, 1024` (cap set by
``recovery_max_substeps``, 4 to 65536). Inside the interval, prescribed position, velocity, and
acceleration follow one quintic (:math:`C^2`) trajectory that lands exactly on the target
values at :math:`t + \Delta t`; host fluid fields and the start-of-interval added mass are held.
Every attempt starts from the complete committed state, including constitutive history. The
host time grid and output times are unchanged, and the run summary reports how many intervals
were subdivided and the largest subdivision used. If the finest subdivision fails, the
interval-start state is restored and the run stops with a named error. A rigid translation of
the line, however fast, changes no chord and no tangent, so it never triggers a subdivision.

``EI = 0`` lines stepped together as one line system (lines joined at ``Free``/``Connect``
points or objects, and the OpenFAST and C-API coupled routes) recover a step that stalls by
halving instead: the interval is re-advanced as two half steps, recursively down to
:math:`\Delta t/64`, with the prescribed motion at each midpoint interpolated linearly between
the interval-start and target position, velocity, and acceleration. If the halves do not
complete the interval, the interval-start state is restored and the full step is re-run, so no
partly subdivided state is ever committed.

**Step size.** Being implicit, the scheme is not bound by the explicit axial-wave limit. In the
snap-load verification a slack line snapping to about 10 MN integrates stably at
:math:`\Delta t = 20\,\mathrm{ms}`, well above the explicit limit of 3.2 ms, with the peak
within 2.7% of a 2 ms reference (:doc:`validation`). The time step must still resolve the
response of interest and, for explicitly integrated points and bodies (below), their own
stability limit.

Attached objects and line ends
------------------------------

Deck ``POINTS``, ``BODIES``, and ``RODS`` records follow the MoorDyn object model
(:ref:`Hall, 2020 <ref-hall2020>`) and are translated into boundary objects of the line system.

- **Fixed, coupled, and vessel points** prescribe the position, velocity, and acceleration of
  line ends.
- **Free points** (``Point3`` buoys, clump weights, connectors) have three translational DOFs.
  Each carries its weight and buoyancy :math:`(\rho V - m)\,g\,\mathbf{e}_z`, point loads, and,
  while below the surface, lumped Morison terms
  :math:`\tfrac{1}{2}\rho\,C_dA\,\lvert\mathbf{u}_r\rvert\mathbf{u}_r + \rho V(1 + C_a)\,\dot{\mathbf{u}}_f`
  with added mass :math:`\rho V C_a`. Points are advanced around the implicit line step with
  implicit-Euler kinematics, and the attached line ends' own blocks are linearly implicit: the
  structural mass diagonal, the Morison added-mass block, and the end-segment elastic tangent
  and axial damping. A light or massless connector is therefore stable at any
  :math:`\Delta t`, including under added mass that exceeds the structural mass, and the end
  segments impose no :math:`\sqrt{k_\text{line}/m_\text{point}}` limit. The coupling to the
  first interior node, drag, and seabed contact at the end node stay at the step-start
  state; the step fails with a named error if a point's motion becomes non-finite.
- **Rigid6 bodies** have six DOFs, a diagonal inertia about the centre of gravity, and the
  mass matrix about the reference point with the centre-of-gravity offset :math:`\mathbf{c}`
  (coupling blocks :math:`\pm m[\mathbf{c}]_\times` and :math:`-m[\mathbf{c}]_\times^2`, and the
  terms :math:`-m\,\boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{c})` and
  :math:`m\,\boldsymbol{\omega}\times(\mathbf{c}\times(\mathbf{c}\times\boldsymbol{\omega}))`).
  The weight acts at the centre of gravity and the buoyancy :math:`\rho g \phi V` at the centre
  of buoyancy. The submerged fraction :math:`\phi` is that of an equivalent sphere of volume
  :math:`V` at the centre of buoyancy, :math:`r_e = (3V/4\pi)^{1/3}`, immersed to
  :math:`h = \min(2r_e, \max(0, \eta - z_b + r_e))`:
  :math:`\phi = h^2(3r_e - h)/(4r_e^3)`, zero when dry and continuous through the surface
  (``bodyWetting moordyn``: :math:`\phi = 1`, as in MoorDyn). A body with linear hydrostatic
  restoring :math:`C_{33}`, :math:`C_{44}`, :math:`C_{55}` (a surface-piercing hull) keeps
  :math:`\phi = 1`; its restoring moment is
  :math:`-C_{44}\tau_1\mathbf{h}_1 - C_{55}\tau_2\mathbf{h}_2`, with :math:`\boldsymbol{\tau}` the
  tilt rotation vector of the body axis that was vertical at the reference pose,
  :math:`\mathbf{h}_1` the horizontal direction of the body x axis and
  :math:`\mathbf{h}_2 = \mathbf{e}_z\times\mathbf{h}_1`, so it does not depend on the heading.
  Lumped Morison drag on the relative velocity (also in still water), fluid inertia (omitted with
  ``bodyHydro moordyn``), and isotropic translational added mass :math:`\rho\phi V C_a` act at
  the reference point, each scaled by :math:`\phi`, and seabed contact acts over a fixed 1 m²
  reference area. The line-end forces act at their attachment offsets with lever moments.
  Accelerations follow the Newton–Euler equations about the reference point; the orientation
  is advanced on SO(3).
- **Rods** are rigid cylinders with a six-DOF rigid state, advanced in the same way (or
  prescribed, for coupled and vessel rods). Their loads are integrated over the ``NumSegs``
  segments of the axis. A cross-section normal to the axis :math:`\mathbf{q}` at angle
  :math:`\varphi` to the vertical is wet below the local free surface over the circular segment
  :math:`A_w(h_w) = r^2\arccos(-h_w/r) + h_w\sqrt{r^2 - h_w^2}`,
  :math:`h_w = (\eta - z)/\sin\varphi`, whose centroid lies
  :math:`-\tfrac{2}{3}(r^2 - h_w^2)^{3/2}/A_w` from the section centre along
  the in-section up-slope direction. With :math:`\eta` linear between the segment end stations,
  each segment splits into fully wet, partially wet, and dry parts at known positions; the wet
  part is integrated with two Gauss points and the partially wet part with ten Gauss points in
  :math:`\theta`, :math:`h_w = r\sin\theta`. The buoyancy :math:`\rho g A_w` per unit length acts
  vertically at the wet centroid, so the hydrostatic force and moment are those of the displaced
  volume of the cylinder, end caps included, at any tilt: a surface-piercing rod carries the exact
  waterplane restoring (second moment :math:`\pi d^4/64`), and the loads are continuous in the
  submergence. Each station also carries, with its wet fraction :math:`f = A_w/A`, the drag
  :math:`\tfrac{1}{2}\rho d f [C_d\lvert\mathbf{u}_n\rvert\mathbf{u}_n + \pi C_{d,ax}
  \lvert\mathbf{u}_t\rvert\mathbf{u}_t]`, the fluid inertia, and the directional added mass
  :math:`\rho A f [C_a(\mathbf{I} - \mathbf{q}\mathbf{q}^T) + C_{a,ax}\mathbf{q}\mathbf{q}^T]`,
  assembled at its offset into a 6×6 added-mass matrix about the rod centre (rotational block
  :math:`C_a \rho A \int s^2\,\mathrm{d}s` about the transverse axes). Each wet end carries the
  MoorDyn end terms: added mass :math:`\rho C_{a,end} V_\text{end}\mathbf{q}\mathbf{q}^T`
  (:math:`V_\text{end} = \tfrac{2}{3}\pi r^3`), axial drag
  :math:`\tfrac{1}{2}\rho C_{d,end} A\lvert u_q\rvert u_q\mathbf{q}`, axial fluid inertia, and,
  in deck waves, the linear dynamic pressure :math:`-(p_\text{dyn} - \rho g\eta^*)A\mathbf{n}_e`
  on the end cap (:math:`\eta^*` the waterline elevation of a surface-piercing rod, zero for a
  submerged one), which carries the axial Froude–Krylov force. Seabed contact is distributed over
  :math:`\max(20, \text{NumSegs}) + 1` stations with :math:`k_\text{Bot}\,d\,\Delta l`.
- **Rods fixed to a body** move with it: with :math:`\mathbf{r}` from the body reference point to
  the rod centre and :math:`\mathbf{T} = [[\mathbf{I}, -[\mathbf{r}]_\times], [\mathbf{0}, \mathbf{I}]]`,
  the rod's structural and added mass enter the body mass matrix as
  :math:`\mathbf{T}^T\mathbf{M}_c\mathbf{T}`, its loads as :math:`(\mathbf{F},
  \mathbf{M} + \mathbf{r}\times\mathbf{F})`, with the velocity terms
  :math:`-\mathbf{M}_c[\boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{r}); \mathbf{0}]`
  and the rod's gyroscopic moment. A **pinned rod** turns about its fixed End A: with the arm
  :math:`\mathbf{r}` from the pin to the centre its acceleration is
  :math:`\mathbf{a} = \boldsymbol{\alpha}\times\mathbf{r} + \boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{r})`,
  and the 6×6 equations are reduced to the three rotational unknowns.
- A **rod pinned to a body** keeps its End A on the body's pin point :math:`\mathbf{p}` (body
  frame) and turns with its own rotation. In reduced coordinates the unknowns are the body twist
  acceleration :math:`\mathbf{x}_b` and the rod's angular acceleration :math:`\mathbf{y}`; the rod's
  centre twist acceleration is
  :math:`\mathbf{X} = \mathbf{G}_b\mathbf{x}_b + \mathbf{G}_c\mathbf{y} + \mathbf{c}` with
  :math:`\mathbf{G}_b = [[\mathbf{I}, -[\mathbf{R}\mathbf{p}]_\times], [\mathbf{0}, \mathbf{0}]]`,
  :math:`\mathbf{G}_c = [[-[\mathbf{l}]_\times], [\mathbf{I}]]` (:math:`\mathbf{l}` from the pin to
  the rod centre) and
  :math:`\mathbf{c} = [\boldsymbol{\omega}_b\times(\boldsymbol{\omega}_b\times\mathbf{R}\mathbf{p})
  + \boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{l}); \mathbf{0}]`. The equations are
  the virtual work of body and rods,
  :math:`\mathbf{M}_b\mathbf{x}_b - \mathbf{F}_b + \sum_k\mathbf{G}_{b,k}^T(\mathbf{M}_k\mathbf{X}_k - \mathbf{F}_k) = \mathbf{0}`
  and :math:`\mathbf{G}_{c,k}^T(\mathbf{M}_k\mathbf{X}_k - \mathbf{F}_k) = \mathbf{0}` for each pinned rod
  :math:`k`: the pin passes the rod's force to the body and no moment, and the rod's rotation
  feels the pin acceleration. (MoorDyn adds every body-attached rod, pinned or not, to the body's
  equations with its full force, moment and mass, and its pinned-rod rotation does not include
  the pin acceleration, so the two formulations can differ for a pinned rod on an accelerating
  body.)
- A **zero-length rod** (``NumSegs`` 0) is a point-like connector: it has no mass, volume or side
  loads, and its end terms act along an axis it does not have, so it reduces to a massless point
  at End A (``Free``, ``Fixed``, ``Coupled`` or on a body).

**Rigid-body initial condition.** Free Rigid6 bodies, free rods and pinned rods start at their
static equilibrium, solved jointly with the ``Free``/``Connect`` points. The unknowns are each
body's position and a rotation-vector increment on SO(3) (for a free rod only the two rotations
normal to its axis), each pinned rod's two rotations about its pin, and each point's position.
A pinned rod's residual rows are its moment about the pin; its net force acts on its parent body
at the pin (a rod pinned to a fixed point has none). A deck with pinned rods takes the
central-difference Jacobian of the full residual. A finite-EI line enters with its axial
stiffness and weight only. The residual is the net force and moment of the loads the
dynamic step applies at rest: weight, buoyancy and hydrostatic restoring, steady-current drag,
and seabed contact, plus the end force of every attached line, re-solved to its static
equilibrium for the moved ends. Moments are divided by the object's length scale (its largest
attachment arm or its half length) and rotations multiplied by it, so forces and moments carry
comparable weight. A trust-region Newton iteration drives the residual below :math:`10^{-9}` of
the largest weight or line force; the first trust radius is the first Newton step. Its Jacobian
takes the objects' own loads by central differences and each line's contribution analytically,
from the line's static tangent condensed onto its two end nodes,
:math:`\mathbf{K}_c = \mathbf{K}_{ee} - \mathbf{K}_{ei}\mathbf{K}_{ii}^{-1}\mathbf{K}_{ie}`
(axial and geometric stiffness of the tension-only law plus the seabed penalty). Only the lines
attached to a moving object are re-solved. Each re-solve starts from the line's previous shape
rotated and stretched onto the new chord, and falls back to the cold static solve if that start
does not converge. At a slack line whose condensed stiffness vanishes, the iteration switches to
a central-difference Jacobian of the full residual.

A step is accepted when it lowers the residual, or when the Newton step from the trial point is
at most half the current one. The second test lets an object swing along a curved valley: a
buoy on a stiff taut line moves on a sphere about the anchor, and every straight step
overstretches the line by its sagitta, so the force residual alone would reject any step
longer than a few centimetres. After each accepted step the rotation unknowns are folded into
the object's orientation and reset to zero, so the analytic Jacobian stays exact however far
the object has turned. A single step turns an object by at most 0.25 rad. A converged
equilibrium must be stable: the symmetric part of the stiffness
:math:`-\partial\mathbf{R}/\partial\mathbf{x}` must be positive semi-definite. If the solve fails or ends on an unstable
equilibrium (for example a body upside down or yawed against its lines), it restarts from the
deck pose with the rotations held until the objects balance in force, then releases them. The
solve stops with an error naming the object if the residual stays above :math:`10^{-6}` of that
force or the equilibrium it reaches is unstable. The
hydrostatic restoring keeps the deck pose as its reference. ``deck bodyIC`` starts from the deck
pose instead, and bodies prescribed by a ``motionFile`` keep their prescribed pose.

**Rigid-body time integration.** Free Rigid6 bodies, Point3 buoys, free and pinned rods and
``Free``/``Connect`` points are stepped monolithically with their ``EI = 0`` lines
(``bodyScheme monolithic``, the default; a deck without lines only when it names it): one implicit
generalised-alpha step in which the objects take every load, the line reactions among them, at
the lines' force level :math:`t_{n+1-\alpha_f}`, :math:`\alpha_f = \rho_\infty/(1 + \rho_\infty)`.
The unknowns are the end-of-step accelerations :math:`\mathbf{x}` of every such junction (a
body's :math:`(\mathbf{a}, \boldsymbol{\alpha})`, a pinned rod's angular acceleration about its
pin, a point's :math:`\mathbf{a}`, a Point3 buoy being the point it is attached at). The bodies'
and rods' inertia sits at the same level (:math:`\alpha_m = \alpha_f`), so their positions and
velocities follow Newmark with :math:`\beta = 1/4`, :math:`\gamma = 1/2`: they have the principal
roots of the trapezoidal rule, second order and free of numerical dissipation (an undamped
linear body conserves its energy exactly), while the lines keep the dissipation of their
:math:`\rho_\infty`; the lines' own :math:`\alpha_m` would add a period error of order
:math:`(1 - \rho_\infty)(\omega\Delta t)^2` to the body modes. The line-end nodes a body or rod
drives follow its Newmark and inertia level inside the line step as well, so their inertia
enters the junction equation consistently. ``Free``/``Connect`` points and Point3 buoys, whose
inertia is largely that of the line ends and whose stiff modes are not sub-stepped, keep the
lines' own generalised-alpha parameters and dissipation. Rotations are on SO(3):
:math:`\Delta\boldsymbol{\theta} = \Delta t\,\boldsymbol{\omega}_n + \Delta t^2[(\tfrac{1}{2} - \beta)\boldsymbol{\alpha}_n + \beta\boldsymbol{\alpha}_{n+1}]`,
:math:`\mathbf{R}_{n+1} = \exp([\Delta\boldsymbol{\theta}]_\times)\mathbf{R}_n`. Each junction's
equation of motion is evaluated at the generalised-alpha level: positions, velocities and
accelerations at :math:`t_{n+1-\alpha_f}` (the orientation on the geodesic between
:math:`\mathbf{R}_n` and :math:`\mathbf{R}_{n+1}`), with the object's mass and
added mass, gyroscopic terms, hydrostatics, drag and seabed contact, and the line loads.

The iteration is a multilevel Newton. For a trial :math:`\mathbf{x}` every attached line is
stepped from its committed state with the resulting end motion prescribed, and returns its
reaction at the alpha level, :math:`\mathbf{f}_e = -\mathbf{R}_E` (the end rows of its own
residual :math:`\mathbf{M}\mathbf{a}_{\alpha_m} + \mathbf{f}_\text{int}(\mathbf{q}_{\alpha_f}) - \mathbf{f}_\text{ext}`),
so the end-node inertia, drag and added mass enter exactly, and its condensed dynamic end
stiffness :math:`\mathbf{K}_e = \mathbf{S}_{EE} - \mathbf{S}_{EF}\mathbf{S}_{FF}^{-1}\mathbf{S}_{FE}`
of the step's effective tangent :math:`\mathbf{S}`, the exact implicit-function derivative of
the reaction with the interior re-solved. The junction Jacobian is the objects' own part
(finite differences of each object's residual, which carry the mass, gyroscopic, hydrodynamic
and contact terms and the turning lever arms) plus
:math:`\beta\Delta t^2\,\boldsymbol{\Gamma}^T\mathbf{K}\boldsymbol{\Gamma}` from the lines, with
:math:`\boldsymbol{\Gamma}` the map from the junction accelerations to the attachment positions. A
Newton step turns an object by at most 0.25 rad. The line re-solves of later iterations start
from the previous solve shifted by the linear map :math:`\mathbf{S}_{FF}^{-1}\mathbf{S}_{FE}` and
reuse its factorisation; the first solve of a step starts from the generalised-alpha predictor
corrected the same way. The condensed stiffness is refreshed every twenty steps, after a step that
needed more than two iterations, and from the third iteration on. The step converges when every
junction residual is within :math:`10^{-6}` of its object's force scale (moments times its size);
the predictor itself is accepted only within :math:`10^{-9}`, so a body at rest does not drift.
Typical steps take two iterations. A step that does not converge in twelve is retried as two
half steps, to a depth of six, and then stops naming the object with the largest residual.

The scheme is implicit in every term, so it needs no sub-stepping for stability: a moored buoy
and a rod spar released out of balance run at :math:`\Delta t` = 0.1 s and 0.5 s. Accuracy is set by
:math:`\Delta t` against the modes of interest: the object period error is the trapezoidal
:math:`(\omega\Delta t)^2/12`, the smallest error constant of any unconditionally stable
second-order method (Dahlquist's second barrier), and twice the :math:`(\omega\Delta t)^2/24` of the
explicit central difference of the staggered scheme. By default (``accuracy bodySubstep``) each
coupling step is therefore divided into equal sub-steps with
:math:`\omega\,\Delta t_\text{sub} \le 0.28 \approx 0.4/\sqrt{2}` for the stiffest body-mooring
frequency estimate below (bodies and rods), which matches the phase accuracy of the staggered
scheme's own sub-steps (:math:`\omega\,\Delta t_\text{sub} \le 0.4`); ``none bodySubstep`` takes the
coupling step as given. The verification tests (CTest ``monolithic_gates``) check: energy of an undamped
linear spring-mass body conserved to :math:`5\times10^{-10}` at :math:`\Delta t` = 0.1 s for any
:math:`\rho_\infty`; observed order 2 in positions and fairlead tensions on the moored buoy and in
surge, heave and pitch on the V-SP spar; gyroscopic precession of a torque-free symmetric top
at second order; and an error no larger than the staggered scheme's at :math:`\Delta t` = 0.1 s and
0.05 s on the moored buoy and rod spar. A finite-EI cable pinned to its object is stepped once
per step with the attachment motion of the Newmark predictor, and its end force is held
through the iterations; a cable clamped to its object is re-stepped at every outer iteration
(*Cables clamped to bodies and rods*, below). Standalone decks with a ``motionFile``, a
``FAILURE`` section, or Coupled/Vessel rods, and the OpenFAST and C-API coupled routes, use the
staggered scheme.

``bodyScheme staggered`` selects the staggered central-difference predictor-corrector
(Newmark :math:`\gamma = 1/2`, :math:`\beta = 0`). The body is first moved to
:math:`\mathbf{r}_{n+1} = \mathbf{r}_n + \Delta t\,\mathbf{v}_n + \tfrac{1}{2}\Delta t^2\mathbf{a}_n`
(the rotation likewise on SO(3)). The lines are then stepped to :math:`t_{n+1}` with that end
motion. The end-of-step acceleration :math:`\mathbf{a}_{n+1}` solves the 6-DOF equation of
motion with the resulting line loads and every body load at :math:`t_{n+1}`, and the velocity
is the trapezoidal average :math:`\mathbf{v}_n + \tfrac{1}{2}\Delta t(\mathbf{a}_n + \mathbf{a}_{n+1})`.
Two terms of that solve are implicit, both exact:

- the inertia of the attached line-end nodes, which keeps a body lighter than the line
  inertia it carries free of the partitioned added-mass instability;
- the velocity derivative of the body drag and seabed damping.

The line loads therefore include the lines' own response over the step (waves and drag on the
lines) without a one-step lag. The scheme is second order and adds no numerical damping, so a
moored body's response converges with :math:`\Delta t`. The line position and velocity
derivatives at the attachment are not used: they describe the end segment with its neighbouring
node held fixed, and so overstate the stiffness and damping the body actually feels. The position
update is explicit. Each coupled step is therefore divided into equal sub-steps with
:math:`\omega\,\Delta t_\text{sub} \le 0.4`, where :math:`\omega` estimates the stiffest
body-mooring frequency. Each attached line acts as the directional spring
:math:`k_a\mathbf{t}\mathbf{t}^T + (T/L)(\mathbf{I} - \mathbf{t}\mathbf{t}^T)` (:math:`\mathbf{t}`
the chord direction, :math:`L` the chord length, :math:`k_a = 1/\sum l_0/EA`, :math:`T` the end
tension), carried to the reference point and added to the hydrostatic restoring; the inertia is
the structural plus added mass. The estimate omits the line-end inertia and the coupling between
the diagonal terms, so it can understate the coupled frequency by about a factor 1.6; stability
alone (the central-difference bound :math:`\omega\,\Delta t < 2`) would allow a limit near 1.
The limit of 0.4 is set by accuracy: a rod released out of balance at :math:`\Delta t = 0.1` s
follows the fairlead-tension transient of a :math:`\Delta t = 0.005` s run within 1.5 %.
Both schemes start from the instantaneous acceleration of the initial state.

**Mixed topologies.** One step marches every object of a deck together: bodies with their fixed
and pinned rods, free and pinned rods, the ``EI = 0`` lines with their ``Free``/``Connect``
points, and the finite-EI cables. In the monolithic scheme the points are junctions of the
Newton iteration above. In the staggered scheme they are integrated on their summed end loads
within the line step. A cable's End A is driven with its attachment's predicted motion and its
End B held; its end force enters the object's equation as a point load at the attachment
(pinned: force only). The cable end-node inertia is not implicit, unlike the ``EI = 0`` line
ends.

**Cables clamped to bodies and rods.** A ``Rigid`` or spring end connection at a cable's End A
on a body point or a rod end fixes the connection direction in the object,
:math:`\mathbf{d} = \mathbf{R}_j\mathbf{d}_0`, so the end tangent follows the object's
rotation; only the direction of :math:`\mathbf{m}` is prescribed and its magnitude (the axial
stretch) stays free. The cable returns to the object its end force at the attachment and the
connection moment :math:`\mathbf{M}_e = \mathbf{m}\times\mathbf{g}_m`, with
:math:`\mathbf{g}_m` the constraint (or spring) force on the tangent coordinates. In the
monolithic step the cable is re-stepped from the step entry at every outer iteration with the
iterate's end-of-step motion, orientation, angular velocity and angular acceleration, and its
end force and moment load the object at the force level
:math:`(1-\alpha_f)\,\mathbf{f}_{n+1} + \alpha_f\,\mathbf{f}_n`. Its condensed dynamic end
stiffness with respect to the object's translation and rotation enters the Newton matrix by
differences of the cable step. On a rod end the direction must lie along the rod axis, so the
moment has no component about the axis, which carries no degree of freedom. The static solve
carries the cable first without bending; passes then add the difference between the cubic-Hermite
cable's end force and that surrogate, and the Hermite connection moment, at the previous
pass's pose, until these loads change by less than :math:`10^{-7}` of the force scale. At the
fixed point the object balances the Hermite cable exactly.

**End connections (Hermite lines).** A line end may carry an isotropic rotational spring with
potential :math:`U = \tfrac{1}{2}k\,\theta^2`, where :math:`\theta` is the angle between the end
tangent :math:`\mathbf{m}` and the reference direction :math:`\mathbf{d}_0`. The spring restrains
only the direction of :math:`\mathbf{m}`; its magnitude, which carries the axial stretch, is
free. :math:`k = 0` is a pinned end. A ``Rigid`` end is an exact two-coordinate constraint
that keeps :math:`\mathbf{m}` on the ray of :math:`\mathbf{d}_0`, not a large stiffness. At a coupled
end, :math:`\mathbf{d}_0` rotates with the supporting body and the connection moment is returned
to it. The connection moment equals the bending moment of the line at its end node.

**Prescribed vessel motion.** In the standalone driver a rigid vessel can carry the coupled
points. With reference point :math:`\mathbf{r}(t)` (``vesselRef`` at the reference pose) and
rotation :math:`\mathbf{R}(t)` from vessel to global axes, a point at the vessel-frame offset
:math:`\mathbf{p}` has

.. math::

   \mathbf{x} = \mathbf{r} + \mathbf{R}\mathbf{p}, \qquad
   \dot{\mathbf{x}} = \dot{\mathbf{r}} + \boldsymbol{\omega}\times\mathbf{R}\mathbf{p}, \qquad
   \ddot{\mathbf{x}} = \ddot{\mathbf{r}} + \boldsymbol{\alpha}\times\mathbf{R}\mathbf{p}
     + \boldsymbol{\omega}\times(\boldsymbol{\omega}\times\mathbf{R}\mathbf{p}),

with :math:`\boldsymbol{\omega}` and :math:`\boldsymbol{\alpha}` the angular velocity and
acceleration in global axes. These are the prescribed position, velocity and acceleration of the
line ends at the point. An end connection at such a point uses
:math:`\mathbf{d}_0(t) = \mathbf{R}(t)\,\mathbf{d}_0` and, for a ``Rigid`` end, the exact rates
:math:`\dot{\mathbf{d}}_0 = \boldsymbol{\omega}\times\mathbf{d}_0` and
:math:`\ddot{\mathbf{d}}_0 = \boldsymbol{\alpha}\times\mathbf{d}_0 + \boldsymbol{\omega}\times\dot{\mathbf{d}}_0`
in the prescribed tangent direction; the tangent magnitude stays a free axial degree of freedom.
Euler angles follow the OrcaFlex vessel convention
:math:`\mathbf{R} = \mathbf{R}_z(\psi)\mathbf{R}_y(\theta)\mathbf{R}_x(\phi)` (roll
:math:`\phi`, pitch :math:`\theta`, yaw :math:`\psi`), for which

.. math::

   \boldsymbol{\omega} = \dot\psi\,\mathbf{e}_z + \dot\theta\,\mathbf{e}_{y'} + \dot\phi\,\mathbf{e}_{x''},
   \qquad \mathbf{e}_{y'} = \mathbf{R}_z\mathbf{e}_y, \quad
   \mathbf{e}_{x''} = \mathbf{R}_z\mathbf{R}_y\mathbf{e}_x,

and :math:`\boldsymbol{\alpha}` is its exact time derivative.

With a displacement RAO table the vessel motion is the linear response to the deck sea. A wave
component of amplitude :math:`a_i` and frequency :math:`\omega_i` has the elevation
:math:`a_i\cos(\omega_i t - \varepsilon_i)` at the reference point, where
:math:`\varepsilon_i = k_i(x_r\cos\beta + y_r\sin\beta) + \varphi_i` carries the component's
phase :math:`\varphi_i` and the reference-point position. Degree of freedom :math:`j` then moves as

.. math::

   \xi_j(t) = r(t)\sum_i a_i A_j(\omega_i, \beta)\cos\big(\omega_i t - \varepsilon_i - P_j(\omega_i, \beta)\big),

with RAO amplitude :math:`A_j` and phase lag :math:`P_j` relative to the wave crest at the
reference point, the relative heading :math:`\beta` equal to the component's own direction of
travel (so each component of a spread or multi-train sea takes the RAO at its heading), and the
``rampTime`` factor :math:`r(t)`. The complex RAO :math:`A_j e^{-iP_j}` is interpolated linearly
in period and heading inside the table; headings are taken modulo 360°, but a heading beyond the
last tabulated one (between 345° and 360° in a 0–345° table, for example) is rejected rather than
interpolated across the wrap to the first. The velocities and accelerations are the exact time
derivatives of
:math:`\xi_j`, and for an irregular sea the motion spectrum is
:math:`S_{\xi_j}(\omega) = A_j^2(\omega)\,S_\eta(\omega)` component by component.

**Discrete attachments (Hermite lines).** A buoyancy module or clump of mass :math:`m`,
displaced volume :math:`V`, normal and axial drag areas :math:`C_dA` and :math:`C_dA_x`, and
added-mass coefficient :math:`C_a` is
lumped at a node :math:`j` of the cable. It adds the constant mass :math:`m + C_a\rho V` to the
three translational DOFs of the node and the nodal load

.. math::

   \mathbf{f}_j = -(m - \rho V)\,g\,\mathbf{e}_z
     + \rho V (1 + C_a)\,\dot{\mathbf{u}}
     + \tfrac{1}{2}\rho\,\big(C_dA\,\lvert\mathbf{w}_n\rvert\mathbf{w}_n
       + C_dA_x\,\lvert\mathbf{w}_t\rvert\mathbf{w}_t\big),

with the fluid velocity :math:`\mathbf{u}` and acceleration :math:`\dot{\mathbf{u}}` at the node and
the relative velocity :math:`\mathbf{w} = \mathbf{u} - \dot{\mathbf{r}}_j` split normal and along
the nodal tangent :math:`\mathbf{m}_j/\lvert\mathbf{m}_j\rvert` (the line's Morison law with
:math:`C_dA = C_{d,n} d` and :math:`C_dA_x = \pi C_{d,t} d`). The drag enters the velocity tangent
and, through :math:`\mathbf{m}_j`, the position tangent exactly. The statics start from the
equilibrium with the net weight smeared over the node's two elements and continue to the lumped
load in warm-started Newton solves; in a deck current the drag at rest is added and the final
solve is repeated with the drag of its own solution until the nodes stop moving. As the pitch :math:`p` of identical modules shrinks, the discrete equilibrium
converges to the smeared section at :math:`O(p^2)`.

**Line failure.** A ``FAILURE`` row detaches a line end from its point at a prescribed time or
when the end tension reaches a threshold, whichever occurs first. The end is transferred onto a
pre-allocated reserve free point that then moves with the attached line ends. On a Rigid6 deck
the bodies keep their staggered step with the remaining lines, and the reserve point is
integrated on its line-end loads within the same line step.

**Active line control.** A ``CONTROL`` row assigns a control channel to an ``EI = 0`` line. The
commanded :math:`\Delta L` and :math:`\Delta\dot L` change the unstretched length and its rate
of the End-A element, :math:`L_0 = L_{0,\text{base}} + \Delta L`; the structural mass and the
point mass shares are updated consistently, while ``BA`` keeps its initial resolution.
Commands that would make the element length non-positive are rejected.

Tension and coupled loads
-------------------------

An interior node tension is the segment tension: the average of the two adjacent element
tensions (``EI = 0``), or the length-weighted average of the element-mean axial forces of the
two adjacent elements (Hermite), including axial damping and constitutive-state contributions.
``FairTen``, ``AnchTen`` and the end-node tension channels report the line-end force instead
(:doc:`conventions`, :doc:`outputs`). By the convention in :doc:`conventions` the reported
tension is the effective tension. The
load returned to a coupled object is the reaction of the cable at the coupled node,
:math:`-(\mathbf{M}\mathbf{a} + \mathbf{f}_\text{int} - \mathbf{f}_\text{ext})` at the prescribed
DOFs, so it includes the line's inertia and hydrodynamic reaction there. Its derivatives with
respect to coupled position, velocity, and acceleration are analytic, with the free DOFs
statically condensed; the acceleration derivative is
:math:`-(\mathbf{M}_{cc} - \mathbf{M}_{cf}\mathbf{M}_{ff}^{-1}\mathbf{M}_{fc})` including added mass
(:doc:`coupling_boundary`).

Linear algebra
--------------

All global operators are stored and factorised in LAPACK band form (``DGBTRF``/``DGBSV``), with
Dirichlet conditions applied in-band and workspaces allocated once, so a dynamic step performs
no heap allocation. The Hermite static energy minimisation also uses the banded Cholesky
factorisation ``DPBTRF``/``DPBTRS``, whose success is the positive-definiteness test of its
shifted tangent. The Hermite operator has half-bandwidth :math:`k_l = k_u = 11` (6 DOFs per
node, 12 per element). The ``EI = 0`` bandwidth is computed from the element connectivity after
removing prescribed DOFs; for a line numbered node by node it is 5. A dense-assembly reference
path packs the assembled matrix into band storage and solves it with ``DGBSV``.
