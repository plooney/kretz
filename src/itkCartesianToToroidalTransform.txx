
#ifndef _itkCartesianToToroidalTransform_txx
#define _itkCartesianToToroidalTransform_txx

#include "itkCartesianToToroidalTransform.h"


namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions>
CartesianToToroidalTransform<TScalarType, NDimensions>::
CartesianToToroidalTransform():Superclass( 0)
{
    return;
}

// Destructor
template<class TScalarType, unsigned int NDimensions>
CartesianToToroidalTransform<TScalarType, NDimensions>::
~CartesianToToroidalTransform()
{
    return;
}

// Print self
template<class TScalarType, unsigned int NDimensions>
void
CartesianToToroidalTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os,indent);
}


template<class TScalarType, unsigned int NDimensions>
double CartesianToToroidalTransform<TScalarType, NDimensions>::Interpolate(double x, std::vector<std::pair<double,double>> table) const {
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

template<class TScalarType, unsigned int NDimensions>
void CartesianToToroidalTransform<TScalarType, NDimensions>::SetTableTheta(const TableType table){
    m_TableAnglesTheta = table;
}

template<class TScalarType, unsigned int NDimensions>
void CartesianToToroidalTransform<TScalarType, NDimensions>::SetTablePhi(const TableType table){
    m_TableAnglesPhi = table;
}

// Transform a point
template<class TScalarType, unsigned int NDimensions>
typename CartesianToToroidalTransform<TScalarType, NDimensions>::OutputPointType
CartesianToToroidalTransform<TScalarType, NDimensions>::
TransformPoint(const InputPointType &point) const 
{

    if (NDimensions != 3) {
        itkExceptionMacro(<< "Method not applicable for dimension not equal to 3.");
        return OutputPointType();
    }

    OutputPointType opoint;
    double x,y,z;

    x = point[0];
    y = point[1];
    z = point[2];

    double pRB, phi, theta;
    double xsq, ysq, zsq, rDsq, tworDz, tmpone, tmptwo, tmpthree, denom;


    phi = std::atan(y/(m_BModeRadius-z));
    double arg1 = std::pow( m_BModeRadius - (y / std::sin(phi)), 2);
    double arg2 = std::pow( x, 2);
    pRB = std::sqrt( arg1 + arg2 );
    theta = std::asin(x/pRB); 
    pRB -= m_SweepRadius;


    if(m_TableAnglesTheta.size() > 0 && m_TableAnglesPhi.size() > 0){
        opoint[0] = pRB/m_Resolution;
        opoint[1] = Interpolate(theta,m_TableAnglesTheta);
        opoint[2] = Interpolate(phi,m_TableAnglesPhi);
    } else {
        opoint[0] = pRB;
        opoint[1] = theta;
        opoint[2] = phi;
    }

    return opoint;
}

} // namespace

#endif
