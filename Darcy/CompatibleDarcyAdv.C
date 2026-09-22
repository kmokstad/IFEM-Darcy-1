// $Id$
//==============================================================================
//!
//! \file CompatibleDarcyAdv.C
//!
//! \date Sep 21 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand for Darcy advection with a div-compatible basis.
//!
//==============================================================================

#include "CompatibleDarcyAdv.h"
#include "DarcyMaterial.h"

#include "CompatibleOperators.h"
#include "EqualOrderOperators.h"
#include "BlockElmMats.h"
#include "FiniteElement.h"
#include "Function.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Vec3Oper.h"


CompatibleDarcyAdv::CompatibleDarcyAdv (unsigned short int n, int torder)
  : DarcyTransport(n,torder)
{
  npv = n+2;
  Darcy::pp = pp;
  Darcy::cc = cc;
  scalarBasis = n+1;
}


LocalIntegral* CompatibleDarcyAdv::getLocalIntegral (const UintVec& nen,
                                                     size_t, bool neumann) const
{
  const size_t nBlock = 5;
  const size_t nBasis = 5;
  BlockElmMats* result = new BlockElmMats(nBlock,nBasis);

  size_t i = 0;
  result->resize(NMAT,NVEC);
  result->redim(qxqx, nen[i++], 1, 1);
  result->redim(qyqy, nen[i++], 1, 2);
  result->redim(qzqz, nsd == 3 ? nen[i++] : 0, 1, 3);

  result->redim(pp, nen[i], 1, -4);
  result->redim(cc, nen[i], 1, 4);
  result->redimOffDiag(qxp, -1);
  result->redimOffDiag(qyp, -1);
  if (nsd == 3)
    result->redimOffDiag(qzp, -1);
  result->redimOffDiag(cqx, 0);
  result->redimOffDiag(cqy, 0);
  if (nsd == 3)
    result->redimOffDiag(cqz, 0);

  result->finalize();
  result->withLHS = true;
  result->rhsOnly = neumann || !calcMats || m_mode >= SIM::RHS_ONLY;
  return result;
}


bool CompatibleDarcyAdv::initElement (const IntVec& MNPC,
                                      const UintVec& elem_sizes,
                                      const UintVec& basis_sizes,
                                      LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  constexpr size_t nScalars = 2;
  const size_t nLevels = this->getNoSolutions();
  elmInt.vec.resize((nsd+nScalars)*nLevels);

  int ierr = 0;
  size_t j = 0;
  for (size_t level = 0; level < nLevels; ++level)
  {
    size_t ofs = 0;
    IntVec::const_iterator first = MNPC.begin();
    for (size_t i = 0; i < nsd; ++i)
    {
      IntVec::const_iterator last = first + elem_sizes[i];
      ierr += utl::gather(IntVec(first,last),0,1,
                          primsol[level],elmInt.vec[j++],ofs,ofs);
      first = last;
      ofs += basis_sizes[i];
    }

    IntVec MNPC2(first,first+elem_sizes[nsd]);
    for (size_t i = 0; i < nScalars; i++)
      ierr += utl::gather(MNPC2,i,nScalars,
                          primsol[level],elmInt.vec[j++],ofs,ofs);
  }

  if (ierr != 0)
    std::cerr << " *** CompatibleDarcyAdv::initElement: Detected " << ierr
              << " node numbers out of range." << std::endl;

  return ierr == 0;
}


bool CompatibleDarcyAdv::evalIntMx (LocalIntegral& elmInt,
                                    const MxFiniteElement& fe,
                                    const TimeDomain& time, const Vec3& X) const
{
  static constexpr std::array<int,3> fluxMatIdx{ qxqx, qyqy, qzqz };
  static constexpr std::array<int,3> concFluxMatIdx{ cqx, cqy, cqz };

  using CompatibleOps = CompatibleOperators::Weak;
  using ScalarOps = EqualOrderOperators::Weak;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  const Vec3 K = mat->getPermeability(X);
  const double mu = mat->getViscosity();

  if (!elMat.A.empty() && calcMats)
  {
    for (size_t i = 0; i < nsd; ++i)
      ScalarOps::Mass(elMat.A[fluxMatIdx[i]], fe,mu/K[i], i+1);

    CompatibleOps::Gradient(elMat.A,fe,{qxp,qyp,qzp});

    const double D = mat->getDispersivity(X);
    ScalarOps::Laplacian(elMat.A[cc],fe,D,false,scalarBasis);

    const double c = this->concentration(elmInt.vec,fe,0);
    for (size_t i = 0; i < nsd; ++i)
      elMat.A[concFluxMatIdx[i]].outer_product(
          fe.grad(scalarBasis).getColumn(i+1),fe.basis(i+1),true,-c*fe.detJxW);
  }

  Vec3 body(gravity);
  if (bodyforce)
    body += (*bodyforce)(X);
  if (!body.isZero())
  {
    body *= mat->getDensity(this->concentration(elmInt.vec,fe,0));
    CompatibleOps::Source(elMat.b,fe,body,{Fqx,Fqy,Fqz});
  }

  if (source)
    ScalarOps::Source(elMat.b[Fp], fe, (*source)(X), 1, scalarBasis);
  if (sourceC)
    ScalarOps::Source(elMat.b[Fc], fe, (*sourceC)(X), 1, scalarBasis);

  if (bdf.getActualOrder() > 0)
  {
    const double phiDt = mat->getPorosity(X)/time.dt;
    double c = 0.0;
    for (int level = 1; level <= bdf.getOrder(); ++level)
      c -= this->concentration(elmInt.vec, fe, level) * phiDt*bdf[level];
    ScalarOps::Source(elMat.b[Fc], fe, c, 1, scalarBasis);
    if (!elMat.A.empty() && calcMats)
      ScalarOps::Mass(elMat.A[cc], fe, phiDt*bdf[0], scalarBasis);
  }

  return true;
}


bool CompatibleDarcyAdv::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X, const Vec3& normal) const
{
  if (!flux && !vflux && !tflux)
  {
    std::cerr << " *** CompatibleDarcyAdv::evalBou: No fluxes." << std::endl;
    return false;
  }

  EqualOrderOperators::Weak::Source(static_cast<ElmMats&>(elmInt).b[Fp],
                                    fe, -this->getFlux(X,normal),
                                    1, scalarBasis);
  return true;
}


bool CompatibleDarcyAdv::evalSol2 (Vector& s, const Vectors& eV,
                                   const FiniteElement& fe, const Vec3& X) const
{
  s = this->darcyFlux(eV,fe).vec(nsd);
  s.push_back(source ? (*source)(X) : 0.0);
  s.push_back(sourceC ? (*sourceC)(X) : 0.0);
  s.push_back(mat->getPorosity(X));

  const Vec3 dC = this->concentrationGradient(eV,fe,0);
  s.push_back(dC.ptr(),dC.ptr()+nsd);
  const Vec3 K = mat->getPermeability(X);
  s.push_back(K.ptr(),K.ptr()+nsd);
  return true;
}


Vec3 CompatibleDarcyAdv::darcyFlux (const Vectors& vec,
                                    const FiniteElement& fe,
                                    size_t level) const
{
  Vec3 result;
  size_t ivec = level*(nsd+2);
  for (size_t i = 0; i < nsd; ++i, ++ivec)
    result[i] = vec[ivec].dot(fe.basis(i+1));
  return result;
}


double CompatibleDarcyAdv::pressure (const Vectors& vec,
                                     const FiniteElement& fe,
                                     size_t level) const
{
  return fe.basis(scalarBasis).dot(vec[level*(nsd+2)+nsd]);
}


double CompatibleDarcyAdv::concentration (const Vectors& vec,
                                          const FiniteElement& fe,
                                          size_t level) const
{
  return fe.basis(scalarBasis).dot(vec[level*(nsd+2)+nsd+1]);
}


Vec3 CompatibleDarcyAdv::concentrationGradient (const Vectors& vec,
                                                const FiniteElement& fe,
                                                size_t level) const
{
  Vector gradient;
  fe.grad(scalarBasis).multiply(vec[level*(nsd+2)+nsd+1],gradient,true);
  return Vec3(gradient);
}


size_t CompatibleDarcyAdv::getNoFields (int fld) const
{
  return fld > 1 ? 3*nsd+3 : nsd+2;
}


std::string CompatibleDarcyAdv::getField1Name (size_t i, const char* pfx) const
{
  if (i == 11)
    return nsd == 2 ? "q_x&&q_y&&p&&c" : "q_x&&q_y&&q_z&&p&&c";
  if (i >= nsd+2u)
    return "";

  static const char* names[5] = {"q_x", "q_y", "q_z", "p", "c"};
  const size_t index = i < nsd ? i : i-nsd+3;
  return pfx ? pfx + std::string(" ") + names[index] : names[index];
}
