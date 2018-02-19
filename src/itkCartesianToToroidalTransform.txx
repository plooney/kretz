
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
void CartesianToToroidalTransform<TScalarType, NDimensions>::SetTable1(const TableType table){
    m_TableAngles1 = table;
}

template<class TScalarType, unsigned int NDimensions>
void CartesianToToroidalTransform<TScalarType, NDimensions>::SetTable2(const TableType table){
    m_TableAngles2 = table;
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

    double pRB, pVAngle, pBAngle;
    double xsq, ysq, zsq, rDsq, tworDz, tmpone, tmptwo, tmpthree, denom;

    /*rDsq = std::pow(m_BModeRadius,2);
    zsq = std::pow(z,2);
    ysq = std::pow(y,2);
    tworDz = 2*m_BModeRadius*z;
    xsq = std::pow(x,2);
    tmpone = rDsq + ysq + zsq -tworDz;
    tmptwo = xsq + ysq + zsq - tworDz + 2*rDsq;
    tmpthree = std::sqrt(tmpone)*2*m_BModeRadius;
    pRB = std::sqrt(tmptwo+tmpthree);
    pBAngle = std::asin(x/pRB);
    //Sign difference from Gordon's thesis. Consistent with segmentations.
    //Need to check if correct with medical image
    pVAngle = -std::asin(y/(m_BModeRadius-pRB*cos(pBAngle)));
    */

    rDsq = std::pow(m_BModeRadius,2);
    zsq = std::pow(z,2);
    ysq = std::pow(y,2);
    tworDz = 2*m_BModeRadius*z;
    xsq = std::pow(x,2);
    tmpone = rDsq + ysq + zsq -tworDz;
    tmptwo = xsq + ysq + zsq - tworDz + 2*rDsq;
    tmpthree = std::sqrt(tmpone)*2*m_BModeRadius;
    pRB = std::sqrt(tmptwo+tmpthree);
    pBAngle = std::asin(x/pRB);
    //Sign difference from Gordon's thesis. Consistent with segmentations.
    //Need to check if correct with medical image
    denom = std::sqrt(std::pow(pRB,2)-std::pow(x,2));
    pVAngle = std::asin(y/(m_BModeRadius-denom));

    double b = pBAngle;
    double v = pVAngle;

    if(m_TableAngles1.size() > 0 && m_TableAngles2.size() > 0){
        pRB -= m_SweepRadius;
        opoint[0] = pRB/m_Resolution;
        opoint[1] = Interpolate(b,m_TableAngles1);
        opoint[2] = Interpolate(v,m_TableAngles2);
    } else {
        opoint[0] = pRB;
        opoint[1] = b;
        opoint[2] = v;
    }

    return opoint;
}

} // namespace

#endif
