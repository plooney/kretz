
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
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetAngles1(const AngleVectorType table){
    m_Angles1 = table;
}

template<class TScalarType, unsigned int NDimensions>
void ToroidalToCartesianTransform<TScalarType, NDimensions>::SetAngles2(const AngleVectorType table){
    m_Angles2 = table;
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

    double b, v, rB;

    b = m_Angles1[point[1]];
    v = m_Angles2[point[2]];
    rB = m_Resolution*point[0] + m_SweepRadius;


    double sinvAngle=std::sin(v);
    double cosvAngle=std::cos(v);

    double sinbAngle=std::sin(b);
    double cosbAngle=std::cos(b);

    x = rB*sinbAngle;
    y = sinvAngle*(rB*cosbAngle-m_BModeRadius);
    z = -m_BModeRadius*(cosvAngle-1)+rB*cosbAngle*cosvAngle;

    opoint[0] = x;
    opoint[1] = y;
    opoint[2] = z;

    return opoint;
}

} // namespace

#endif
