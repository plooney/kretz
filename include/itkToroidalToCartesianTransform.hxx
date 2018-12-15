/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkToroidalToCartesianTransform_hxx
#define itkToroidalToCartesianTransform_hxx

#include "itkToroidalToCartesianTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions>
ToroidalToCartesianTransform<TScalarType, NDimensions>::
ToroidalToCartesianTransform():Superclass( 0)
{
    return;
}

// Destructor
template<class TScalarType, unsigned int NDimensions>
ToroidalToCartesianTransform<TScalarType, NDimensions>::
~ToroidalToCartesianTransform()
{
    return;
}

// Print self
template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetTableTheta(const TableType table){
    m_TableAnglesTheta = table;
}

template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetTablePhi(const TableType table){
    m_TableAnglesPhi = table;
}

template<class TScalarType, unsigned int NDimensions>
double ToroidalToCartesianTransform<TScalarType, NDimensions>::Interpolate(double x, std::vector<std::pair<double,double>> table) const {
  // Assumes that "table" is sorted by .first
  // Check if x is out of bound
  if (x > table.back().first) return std::numeric_limits<double>::max();
  if (x < table[0].first) return -std::numeric_limits<double>::max();
  std::vector<std::pair<double, double> >::iterator it, it2;
  // INFINITY is defined in math.h in the glibc implementation
  it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -std::numeric_limits<double>::max()));
  // Corner case
  if (it == table.begin()) return it->second;
  it2 = it;
  --it2;
  return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}


// Transform a point
template<class TScalarType, unsigned int NDimensions>
typename ToroidalToCartesianTransform<TScalarType, NDimensions>::OutputPointType
ToroidalToCartesianTransform<TScalarType, NDimensions>::
TransformPoint(const InputPointType &point) const 
{
  if (NDimensions != 3) {
      itkExceptionMacro(<< "Method not applicable for dimension not equal to 3.");
      return OutputPointType();
  }

  OutputPointType opoint;
  double x,y,z;

  double theta, phi, rB;

  phi = m_TableAnglesPhi[point[2]].first;
  theta = m_TableAnglesTheta[point[1]].first;
  rB = m_Resolution*point[0] + m_SweepRadius;

  double sinPhi = std::sin(phi);
  double cosPhi = std::cos(phi);

  double sinTheta = std::sin(theta);
  double cosTheta = std::cos(theta);

  x = rB*sinTheta;
  y = sinPhi*(rB*cosTheta -m_BModeRadius);
  z = m_BModeRadius*(1-cosPhi)+rB*cosTheta*cosPhi;

  opoint[0] = x;
  opoint[1] = y;
  opoint[2] = z;

  return opoint;
}

template< typename TScalarType, unsigned int NDimensions >
template< typename TImage >
typename ToroidalToCartesianTransform< TScalarType, NDimensions >::BoundsArrayType
ToroidalToCartesianTransform<TScalarType, NDimensions>
::ComputeBounds(const TImage * image, const Self * t2c)
{
  typedef TImage ImageType;
  MeshType::Pointer mesh = MeshType::New();
  ImageRegionConstIteratorWithIndex<ImageType> imageIterator(image,image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
  {
    typename ImageType::IndexType index = imageIterator.GetIndex();
    typename ImageType::PointType p;
    p[0] = index[0];
    p[1] = index[1];
    p[2] = index[2];
    typename ImageType::PointType p_cart = t2c->TransformPoint(p);

    unsigned long i = mesh->GetNumberOfPoints();
    mesh->SetPoint(i, p_cart);
    ++imageIterator;
  }
  typename MeshType::PointsContainer::Pointer points = mesh->GetPoints();

  typename BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();

  boundingBox->SetPoints(points);
  boundingBox->ComputeBoundingBox();

  BoundsArrayType bounds = boundingBox->GetBounds();

  return bounds;
}

} // namespace itk

#endif
