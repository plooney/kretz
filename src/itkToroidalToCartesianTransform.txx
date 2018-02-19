
#ifndef _itkToroidalToCartesianTransform_txx
#define _itkToroidalToCartesianTransform_txx

#include "itkToroidalToCartesianTransform.h"

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
void
ToroidalToCartesianTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}

template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetTable1(const TableType table){
    m_TableAngles1 = table;
}

template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetTable2(const TableType table){
    m_TableAngles2 = table;
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

    phi = m_TableAngles2[point[2]].first;
    theta = m_TableAngles1[point[1]].first;
    rB = m_Resolution*point[0] + m_SweepRadius;

    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);

    double sinTheta = std::sin(theta);
    double cosTheta = std::cos(theta);

    x = rB*sinTheta;
    y = -sinPhi*(rB*cosTheta -m_BModeRadius);
    z = m_BModeRadius*(1-cosPhi)+rB*cosTheta*cosPhi;

    opoint[0] = x;
    opoint[1] = y;
    opoint[2] = z;

    return opoint;
}

} // namespace

#endif
