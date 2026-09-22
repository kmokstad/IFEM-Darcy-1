// $Id$
//==============================================================================
//!
//! \file CompatibleDarcyAdv.h
//!
//! \date Sep 21 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand for Darcy advection with a div-compatible basis.
//!
//==============================================================================

#ifndef _COMPATIBLE_DARCY_ADV_H_
#define _COMPATIBLE_DARCY_ADV_H_

#include "DarcyTransport.h"


/*!
  \brief Darcy advection integrand using a divergence-compatible flux basis.
  \details The Darcy flux has one scalar basis for each spatial component,
  while pressure and concentration share the final scalar basis.
*/

class CompatibleDarcyAdv : public DarcyTransport
{
  //! \brief Element right-hand-side vector indices.
  enum ResidualVectors { Fqx = 1, Fqy = 2, Fqz = 3, Fp = 4, Fc = 5, NVEC = 6 };

  //! \brief Element tangent matrix indices.
  enum TangentMatrices {
    qxqx =  1, qyqy =  2, qzqz =  3, pp  =  4, cc =  5,
    qxqy =  6, qxqz =  7, qxp  =  8, qxc =  9,
    qyqx = 10, qyqz = 11, qyp  = 12, qyc = 13,
    qzqx = 14, qzqy = 15, qzp  = 16, qzc = 17,
    pqx  = 18, pqy  = 19, pqz  = 20, pc  = 21,
    cqx  = 22, cqy  = 23, cqz  = 24, cp  = 25,
    NMAT = 26
  };

  using IntVec  = const std::vector<int>;    //!< Convenience type alias
  using UintVec = const std::vector<size_t>; //!< Convenience type alias

public:
  //! \brief The constructor initializes the integrand.
  explicit CompatibleDarcyAdv(unsigned short int n, int torder = 0);

  using DarcyBase::getLocalIntegral;
  //! \brief Returns a mixed local integral container.
  LocalIntegral* getLocalIntegral(const UintVec& nen,
                                  size_t, bool neumann) const override;

  using DarcyBase::initElement;
  //! \brief Extracts element vectors from a mixed-basis solution.
  bool initElement(const IntVec& MNPC,
                   const UintVec& elem_sizes, const UintVec& basis_sizes,
                   LocalIntegral& elmInt) override;

  using DarcyBase::evalIntMx;
  //! \brief Evaluates the mixed interior integrand.
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const TimeDomain& time, const Vec3& X) const override;

  using DarcyBase::evalBou;
  //! \brief Evaluates a prescribed normal Darcy flux.
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Avoids the equal-order transport finalization.
  bool finalizeElement(LocalIntegral&) override { return true; }

  //! \brief Evaluates secondary fields at a result point.
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary fields.
  size_t getNoFields(int fld) const override;

  //! \brief Returns the name of a primary field.
  std::string getField1Name(size_t i, const char* prefix) const override;

protected:
  //! \brief Evaluates the Darcy flux at an integration point.
  Vec3 darcyFlux(const Vectors& vec, const FiniteElement& fe,
                 size_t level = 0) const;

  //! \brief Evaluates pressure at an integration point.
  double pressure(const Vectors& vec, const FiniteElement& fe,
                  size_t level) const override;

  //! \brief Evaluates concentration at an integration point.
  double concentration(const Vectors& vec, const FiniteElement& fe,
                       size_t level) const;

  //! \brief Evaluates the concentration gradient at an integration point.
  Vec3 concentrationGradient(const Vectors& vec, const FiniteElement& fe,
                             size_t level) const;

private:
  size_t scalarBasis; //!< Index of basis used for pressure and concentration
};

#endif
