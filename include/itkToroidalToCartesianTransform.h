

#ifndef __itkToroidalToCartesianTransform_h
#define __itkToroidalToCartesianTransform_h


#include <iostream>
#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMatrix.h"
#include <vector>
#include <itkArray.h>
#include <math.h>
#include <limits>


namespace itk
{

/** \brief Toroidal transformation of a vector space (e.g. space coordinates).
 *
 * Transforms first two coordinates form polar space <alpha,radius> to cartesian
 * coordinates. Other dimensions are left unchanges. In fact this is generalized
 * cylindric transform:
 * \f[			x_1 = r cos( \alpha ) \f]
 * \f[			x_2 = r sin( \alpha ) \f]
 * \f[			x_n = x_n, \mbox{ n>=2 } \f]
 *
 *
 * \par
 * Center of the polar transform is a center of coordinate system < 0, 0 >.
 *
 * Dimension must be at least 2 or an exception is thrown during transform.
 *
 * Extent of input in first dimension (alpha) should be only < 0, 2*pi ).
 *
 * \author Pádraig Looney, Nuffield Department of Obstetrics and Gynaecology, University of Oxford.
 *
 * \ingroup Transforms
 */
template <
        class TScalarType=double,          // Data type for scalars (float or double)
        unsigned int NDimensions=3>        // Number of dimensions
class ITK_EXPORT ToroidalToCartesianTransform : public Transform< TScalarType, NDimensions, NDimensions >
{
public:

    /** Standard class typedefs. */
    typedef ToroidalToCartesianTransform Self;
    typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;
    typedef std::vector<std::pair<double, double> > TableType;

    /** Variables specific to geometry*/
    double m_BModeRadius;
    double m_SweepRadius;
    double m_Resolution;

    TableType m_TableAngles1;
    TableType m_TableAngles2;

    /** New macro for creation of through the object factory.*/
    itkNewMacro( Self )

    itkSetMacro(BModeRadius, double)
    itkSetMacro(SweepRadius, double)
    itkSetMacro(Resolution, double)

    /** Run-time type information (and related methods). */
    itkTypeMacro( ToroidalToCartesianTransform, Transform )

    /** Dimension of the domain space. */
    itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
    itkStaticConstMacro(ParametersDimension, unsigned int, 0);

    /** Standard scalar type for this class. */
    typedef typename Superclass::ScalarType ScalarType;

    /** Standard Jacobian container. */
    typedef typename Superclass::JacobianType JacobianType;

    /** Standard vector type for this class. */
    typedef Vector<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    InputVectorType;
    typedef Vector<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    OutputVectorType;

    /** Standard covariant vector type for this class. */
    typedef CovariantVector<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    InputCovariantVectorType;
    typedef CovariantVector<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    OutputCovariantVectorType;

    /** Standard vnl_vector type for this class. */
    typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    InputVnlVectorType;
    typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    OutputVnlVectorType;

    /** Standard coordinate point type for this class. */
    typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    InputPointType;
    typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)>
    OutputPointType;

    typedef typename Superclass::FixedParametersType FixedParametersType;
    typedef typename Superclass::ParametersType ParametersType;

    /** Method to transform a point.
   * This method transforms first two dimensions of a point from polar
   * coordinates <alpha,radius> to cartesian coordinates.
   */
    OutputPointType     TransformPoint(const InputPointType  &point ) const;

    /** Method to transform a vector - not applicable for this type of
      transform. */
    virtual OutputVectorType TransformVector(const InputVectorType &) const
    {
        itkExceptionMacro(<< "Method not applicable for polar transform." );
        return OutputVectorType();
    }

    /** Method to transform a vnl_vector - not applicable for this type of
      transform. */
    virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
    {
        itkExceptionMacro(<< "Method not applicable for polar transform. ");
        return OutputVnlVectorType();
    }

    /** Method to transform a CovariantVector - not applicable for this type of
      transform */
    virtual OutputCovariantVectorType TransformCovariantVector(
            const InputCovariantVectorType &) const
    {
        itkExceptionMacro(<< "Method not applicable for polar transfrom. ");
        return OutputCovariantVectorType();
    }

    /** Compute the Jacobian Matrix of the transformation at one point - not
      applicable for this type of transform */
    virtual void ComputeJacobianWithRespectToPosition( const InputPointType &, JacobianType &)  const
    {
        itkExceptionMacro(<< "Method not applicable for polar transform. ");
        //return this->m_Jacobian;
    }

    /** Compute the Jacobian Matrix of the transformation at one point - not
      applicable for this type of transform */
    virtual void ComputeJacobianWithRespectToParameters( const InputPointType &, JacobianType &)  const
    {
        itkExceptionMacro(<< "Method not applicable for polar transform. ");
        //return this->m_Jacobian;
    }

    virtual void SetTable1(const TableType);
    virtual void SetTable2(const TableType);

    virtual void SetParameters(const ParametersType &){

    }
    virtual void SetFixedParameters(const FixedParametersType &){

    }

protected:
    ToroidalToCartesianTransform();
    ~ToroidalToCartesianTransform();

    /** Print contents of an ToroidalToCartesianTransform. */
    void PrintSelf(std::ostream &os, Indent indent) const;

private:
    ToroidalToCartesianTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    //virtual double Interpolate(double x, TableType table) const;

}; //class ToroidalToCartesianTransform


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkToroidalToCartesianTransform.txx"
#endif

#endif /* __itkToroidalToCartesianTransform_h */
