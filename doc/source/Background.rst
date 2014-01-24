Background
==========

Laplacian-based approaches represent the surface by the so-called differential
coordinates or Laplacian coordinates :cite:`Alexa:2003`, :cite:`SorkineCohenToledo:2003`. These coordinates are obtained
by applying the Laplacian operator to the mesh vertices:

.. math::

  \boldsymbol{ \delta }_i = \Delta_{S}( \boldsymbol{ p }_i ) = - H_i \cdot \boldsymbol{ n_i }

where :math:`H_i` is the mean curvature ( :math:`\kappa_1 + \kappa_2` ) at the
vertex :math:`v_i`.

The deformation can be formulated by minimizing the difference from the input
surface coordinates :math:`\delta_i`. With a continuous formulation, this would
lead to the minimization of the following energy:

.. math::

  \min_{\boldsymbol{p'}} \int_{\Omega} \| \boldsymbol{\Delta p'} - \boldsymbol{\delta} \| du dv

The Euler-Lagrange equation derived:

.. math::

  \Delta^2 \boldsymbol{p'} = \Delta \boldsymbol{\delta}

When considering the input surface as the parameter domain, the Laplace
operator turns out into the Laplace-Beltrami operator :math:`\Delta_S`:

.. math::

  L^2 \boldsymbol{p'} = L \boldsymbol{ \delta }

which can be separated into 3 coordinate components. Then users can add
positional constraints on some vertices:

.. math::

  \boldsymbol{p'}_j = \boldsymbol{c}_j


Note that the positional constraints can either be incorporated as *hard* or
*soft constraints*.

Adding constraints
++++++++++++++++++

Hard constraints formulation
----------------------------

The problem can thus be expressed as follows:

.. math::

  L^2 \boldsymbol{p'} & =  L \boldsymbol{ \delta } \\
  \boldsymbol{p'}_j   & =  \boldsymbol{c}_j

where :math:`\boldsymbol{c}_j` are the positional constraints, or by
constraining the displacement :math:`\boldsymbol{d}_i = \boldsymbol{c}_i -
\boldsymbol{p}_i`.

which can be rewritten as:

.. math::

  \left(
    \begin{matrix}
      L^2 \\
      0 & I_k
    \end{matrix}
  \right) \cdot
  \left(
    \begin{matrix}
      \boldsymbol{d}_1 \\
      \vdots \\
      \boldsymbol{d}_n
    \end{matrix}
  \right)
  =
  \left(
    \begin{matrix}
      \boldsymbol{0} \\
      \vdots \\
      \boldsymbol{0} \\
      \boldsymbol{c}_{n'+1} - \boldsymbol{p}_{n'+1} \\
      \vdots \\
      \boldsymbol{c}_{n} - \boldsymbol{p}_{n}
    \end{matrix}
  \right)

where :math:`L^2 \in \mathbb{R}^{n \times n}` and :math:`I_k` is the :math:`k
\times k` identity matrix. By eliminating rows in the above linear system, we
finally get a :math:`n' \times n'` sparse linear system where the unknown
vector represents the 3D deformation of the :math:`n'` unconstrained vertices.

Soft constraints formulation
----------------------------

In contrast soft constraints correspond to relax the previous constraint by
incorporating the constraints into the enyergy minimization leading to

.. math::

  \min_{\boldsymbol{p'}} \sum_{i=1}^{N} \| \Delta_S( \boldsymbol{p'}_i - \boldsymbol{\delta}_i \|^2 + \lambda \sum_{j=n'+1}^{n} \| \boldsymbol{p'}_j - \boldsymbol{c}_j \|^2

The minimum of this energy is can be found by solving the normal equations:

.. math::

  \left[
    L^t \cdot L + \lambda^2 \cdot
    \left(
      \begin{matrix}
        0 & 0 \\
        0 & I_k
      \end{matrix}
    \right)
  \right]
  \cdot
  \left(
    \begin{matrix}
      \boldsymbol{p}_1 \\
      \vdots \\
      \boldsymbol{p}_n
    \end{matrix}
  \right)
  =
  L^t \cdot
  \left(
    \begin{matrix}
      \boldsymbol{\delta'}_0 \\
      \vdots \\
      \boldsymbol{\delta'}_n
    \end{matrix}
  \right)
  + \lambda^2 \cdot
  \left(
    \begin{matrix}
      \boldsymbol{0} \\
      \vdots \\
      \boldsymbol{0} \\
      \boldsymbol{c}_{n'+1} \\
      \vdots \\
      \boldsymbol{c}_n
    \end{matrix}
  \right)

Depending on :math:`\lambda` solution can be closed to an interpolation of the
constraints (with large values), or an approximation of them (with low values).

Laplacian discretization
++++++++++++++++++++++++

This approach requires a dicretization of the Laplacian operator, and results
would highly depends on it. There exists several variations of the weights used
in the typically used Laplacian discretization. Here we list few of them:

* The uniform weight: :math:`w_i = 1` and :math:`w_{ij} = 1 / N`;
* :math:`w_{i} = 1` and :math:`w_{ij} = \frac{1}{2} ( \cot \alpha_{ij} + \cot \beta_{ij} )`;
* :math:`w_{i} = 1 / A_i` and :math:`w_{ij} = \frac{1}{2} ( \cot \alpha_{ij} + \cot \beta_{ij} )` (where :math:`A_i` is a local area corresponding for :math:`\boldsymbol{v}_i`)

.. bibliography:: refs.bib
