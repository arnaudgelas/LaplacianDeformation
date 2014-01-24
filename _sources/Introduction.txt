Introduction
============

Manipulating and modifying surface while preserving geometric details has been
an active area of research in geometric modeling due to their applications in
design (e.g. computer graphics, animation), but also in biomedical imaging
(e.g. deforming surface based atlas).

Here, we present one surface-based technique, as opposed to free-form
deformations which deform the ambient 3D space, and based on differential
representations (i.e. gradient based representation, Laplacian based
representation, local frame representation). These techniques are becoming more
and more popular over the last past years, most likely due to their robustness,
speed and ease of implementation. The main idea behind these approaches relies
on the use of a representation that focus on local differential properties and
on preserving these ones when deforming. These approaches remain quite
intuitive and preserve local details throughout the deformation. In the rest of
the paper, we will focus on Laplacian based representation and will focus on methods as described in :cite:`BotschSorkine:2008`

