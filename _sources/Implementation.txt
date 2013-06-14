Classes
=======

There are 2 cases inheriting from a single base class (which is not meant to be
instantiated).

Hard Constraints
----------------

For hard constraints, one should first choose the type of sparse solver to be
used. For performance purpose we strongly advise to use
*VNLSparseLUSolverTraits*, but it is possible to use a custom one with external
libraries not delivered in ITK (such as TAUCS, CHOLMOD, etc.).

.. code-block:: c++

  typedef VNLSparseLUSolverTraits< CoordType > SolverType;
  typedef itk::LaplacianDeformationQuadEdgeMeshFilterWithHardConstraints< MeshType, MeshType, SolverType > FilterType;

The next step is then to allocate the filter and set the input mesh.

.. code-block:: c++

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );

Then set the order (note that **we recommend using order 2**, but it is possible to use higher-orders).

.. code-block:: c++

  filter->SetOrder( 2 );

Then set the method to compute the coefficient. Although theoretically it is
recommended to use *itk::ConformalMatrixCoefficients*, it is possible to use
other weights to speed up the deformation, or for other reasons...

.. code-block:: c++

  typedef itk::ConformalMatrixCoefficients< MeshType > CoefficientType;
  CoefficientType coeff;
  filter->SetCoefficientsMethod( &coeff );

Then set the constraints, either by the means of *SetDisplacement*, either by the means of *SetConstrainedNode*

.. code-block:: c++

  MeshType::VectorType null( 0. );
  filter->SetDisplacement( 150, null );

  MeshType::VectorType d( null );
  d[2] = -5;

  filter->SetDisplacement( 2030, d );

Finally call the *Update* method and get the output mesh

.. code-block:: c++

  filter->Update();

  MeshType::Pointer output = filter->GetOutput();


Soft Constraints
----------------

As above (for hard constraints), one should first choose the type of sparse
solver to be used. For performance purpose we strongly advise to use
*VNLSparseLUSolverTraits*, but it is possible to use a custom one with external
libraries not delivered in ITK (such as TAUCS, CHOLMOD, etc.).

.. code-block:: c++

  typedef VNLSparseLUSolverTraits< CoordType > SolverType;
  typedef itk::LaplacianDeformationQuadEdgeMeshFilterWithSoftConstraints< MeshType, MeshType, SolverType > FilterType;

The next step is then to allocate the filter and set the input mesh.

.. code-block:: c++

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );

Then set the order (note that **we recommend using order 1**, but it is
possible to use higher-orders), and the lambda value which balance in between
interpolation and approximation.

.. code-block:: c++

  filter->SetOrder( 1 );
  filter->SetLambda( 1. );

Then set the method to compute the coefficient. Although theoretically it is
recommended to use *itk::ConformalMatrixCoefficients*, it is possible to use
other weights to speed up the deformation, or for other reasons...

.. code-block:: c++

  typedef itk::ConformalMatrixCoefficients< MeshType > CoefficientType;
  CoefficientType coeff;
  filter->SetCoefficientsMethod( &coeff );

Then set the constraints, either by the means of *SetDisplacement*, either by the means of *SetConstrainedNode*

.. code-block:: c++

  MeshType::VectorType null( 0. );
  filter->SetDisplacement( 150, null );

  MeshType::VectorType d( null );
  d[2] = -5;

  filter->SetDisplacement( 2030, d );

Finally call the *Update* method and get the output mesh

.. code-block:: c++

  filter->Update();

  MeshType::Pointer output = filter->GetOutput();


