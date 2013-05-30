#include "itkQuadEdgeMesh.h"
#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"
#include "itkLaplacianDeformationQuadEdgeMeshFilterWithHardConstraints.h"
#include "VNLSparseLUSolverTraits.h"

int main( int argc, char* argv[] )
{
  const unsigned int Dimension = 3;
  typedef double CoordType;
  typedef itk::QuadEdgeMesh< CoordType, Dimension > MeshType;
  typedef itk::VTKPolyDataReader< MeshType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef VNLSparseLUSolverTraits< CoordType > SolverType;
  typedef itk::LaplacianDeformationQuadEdgeMeshFilterWithHardConstraints< MeshType, MeshType, SolverType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetOrder( 2 );

  typedef itk::ConformalMatrixCoefficients< MeshType > CoefficientType;
  CoefficientType coeff;
  filter->SetCoefficientsMethod( &coeff );

  MeshType::VectorType null( 0. );
  filter->SetDisplacement( 150, null );
  filter->SetDisplacement( 2027, null );
  filter->SetDisplacement( 292, null );
  filter->SetDisplacement( 185, null );
  filter->SetDisplacement( 180, null );
  filter->SetDisplacement( 153, null );
  filter->SetDisplacement( 183, null );
  filter->SetDisplacement( 226, null );

  MeshType::VectorType d( null );
  d[2] = -5;

  filter->SetDisplacement( 2030, d );
  filter->SetDisplacement( 1870, d );
  filter->SetDisplacement( 2012, d );

  d[1] = 1;
  filter->SetDisplacement( 1076, d );

  d[1] = -1;
  filter->SetDisplacement( 1077, d );

  d[1] = 0.;
  filter->SetDisplacement( 1062, d );

  filter->Update();

  typedef itk::VTKPolyDataWriter< MeshType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}
