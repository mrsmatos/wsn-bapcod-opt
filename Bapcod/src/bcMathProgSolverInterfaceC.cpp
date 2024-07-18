/**
 *
 * This file bcMathProgSolverInterfaceC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMathProgSolverInterfaceC.hpp"
#include "bcBapcodInit.hpp"
#include "bcCoefAboundC.hpp"
#include "bcControlParameters.hpp"

MathProgSolverInterface::MathProgSolverInterface(BapcodInit* bapcodInit) :
  _ref(0),
  _formCurrentlyLoaded(false),
  _formCurrentlySaved(false),
  _pureLP(true),
  _ncol(0),
  _nrow(0),
  _bapcodInit(bapcodInit)
{
}

MathProgSolverInterface::MathProgSolverInterface(BapcodInit* bapcodInit, const int & refconst,
                                                 const std::string & name):
  _ref(refconst),
  _formCurrentlyLoaded(false),
  _formCurrentlySaved(false),
  _pureLP(true),
  _ncol(0),
  _nrow(0),
  _bapcodInit(bapcodInit)
{
}

MathProgSolverInterface::~MathProgSolverInterface()
{
}

const int & MathProgSolverInterface::ref() const
{
  return _ref;
}

const long & MathProgSolverInterface::ncol() const
{
  return _ncol;
}

const long & MathProgSolverInterface::nrow() const 
{
  return _nrow;
}

BapcodInit & MathProgSolverInterface::bapcodInit() const
{
  return *_bapcodInit;
}

BapcodInit * MathProgSolverInterface::bapcodInitPtr() const
{
  return _bapcodInit;
}

const ControlParameters& MathProgSolverInterface::param() const
{
  return _bapcodInit->param();
}

ControlParameters& MathProgSolverInterface::param()
{
  return _bapcodInit->param();
}

void MathProgSolverInterface::setLazyConstraintsCallback(MasterConf * mastConfPtr)
{
    std::cerr << "Lazy constraint callback cannot be set!" << std::endl;
}

void MathProgSolverInterface::removeLazyConstraintsCallback()
{
}
