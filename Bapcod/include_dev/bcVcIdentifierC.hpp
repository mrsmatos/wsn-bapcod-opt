/**
 *
 * This file bcVcIdentifierC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCVCIDENTIFIERC_HPP_
#define BCVCIDENTIFIERC_HPP_

#ifdef _WINDOWS
#include <stdint.h>
#endif

struct VcId
{
#ifdef _WINDOWS
  enum VcIdentifier : int64_t
#else
  enum VcIdentifier
#endif
    {
	VarConstrMask                            = 0x01,
	AggregateVariableMask                    = 0x02,
	//MasterVarMask                            = 0x04,
	ArtificialVarMask                        = 0x08,
	VariableMask                             = 0x10               | VarConstrMask,
	InstanciatedVarConstrMask                = 0x20,
	MasterConstrMask                         = 0x40,
	ConstraintMask                           = 0x80               | VarConstrMask,
	MastColumnMask                           = 0x100              | AggregateVariableMask //| MasterVarMask
											                      | VariableMask,
	GlobalArtificialVarMask                  = 0x200              | ArtificialVarMask | VariableMask,
	LocalArtificialVarMask                   = 0x400              | ArtificialVarMask | VariableMask,
	InstanciatedVarMask                      = 0x800              | VariableMask | InstanciatedVarConstrMask,
	OvfVarMask                               = 0x1000             | VariableMask,
	InstanciatedConstrMask                   = 0x2000             | InstanciatedVarConstrMask | ConstraintMask,
	BranchingConstrBaseTypeMask              = 0x8000,
	OvfConstrMask                            = 0x10000            | ConstraintMask,
    MissingColumnMask                        = 0x20000            | MastColumnMask | ArtificialVarMask,
	InstanciatedAggregateVarMask             = 0x40000            | AggregateVariableMask | InstanciatedVarMask,
	InstMasterVarMask                        = 0x80000            | InstanciatedVarMask /*| MasterVarMask */,
	SubProbVariableMask                      = 0x100000           | InstanciatedVarMask,
	SpSetupOvfVarMask                        = 0x200000           | OvfVarMask,
	Base4NonLinearConstraintMask             = 0x400000,
	InstMasterConstrMask                     = 0x800000           | MasterConstrMask | InstanciatedConstrMask,
	//SeparationSpConstrMask                   = 0x1000000,
	InstanciatedBranchingConstrMask          = 0x2000000          | InstanciatedConstrMask
											                      | BranchingConstrBaseTypeMask,
	InstSubProbBranchingConstrMask           = 0x4000000          | BranchingConstrBaseTypeMask
											                      | InstanciatedConstrMask,
	SpecificOvfConstrMask                    = 0x8000000          | OvfConstrMask,
	NonLinearInstConstrMask                  = 0x10000000         | InstanciatedConstrMask
											                      | Base4NonLinearConstraintMask,
	NonLinearInstMastConstrMask              = 0x20000000         | Base4NonLinearConstraintMask | InstMasterConstrMask,
	InstMastConvexityConstrMask              = 0x40000000         | InstMasterConstrMask,
	//InstMastOptimalityCutOffValueConstrMask  = 0x80000000         | InstMasterConstrMask,
	//InstMastArtVarNormalizationConstrMask    = 0x100000000        | InstMasterConstrMask,
	InstMasterBranchingConstrMask            = 0x200000000        | InstMasterConstrMask | BranchingConstrBaseTypeMask,
	RyanAndFosterInstSubProbBranchConstrMask = 0x400000000        | InstSubProbBranchingConstrMask,
	SpVarLbOvfConstrMask                     = 0x800000000        | SpecificOvfConstrMask,
	SpVarUbOvfConstrMask                     = 0x1000000000       | SpecificOvfConstrMask,
	SpLbOvfConstrMask                        = 0x2000000000       | SpecificOvfConstrMask,
	CompSetInstMastBranchConstrMask          = 0x4000000000       | Base4NonLinearConstraintMask
											                      | InstMasterBranchingConstrMask,
	GenVarInstMastBranchConstrMask           = 0x8000000000       | InstMasterBranchingConstrMask,
	BasicConstrInstMastBranchingConstrMask   = 0x10000000000      | InstMasterBranchingConstrMask,
	LocalMasterBranchingConstrMask           = 0x20000000000      | MasterConstrMask | ConstraintMask
											                      | BranchingConstrBaseTypeMask,
	SolMinimalMastColMask                    = 0x40000000000      | VariableMask,
	//UncertainConstrMask					     = 0x80000000000,
    
    //Generic var / constr
	GenVarConstrMask                         = 0x100000000000,
	GenericConstrMask                        = 0x200000000000     | GenVarConstrMask,
	DynamicVarConstrMask                     = 0x400000000000,     
	DynamicGenericConstrMask                 = 0x800000000000     | DynamicVarConstrMask | GenericConstrMask,
	GenericBranchingConstrMask               = 0x1000000000000    | DynamicGenericConstrMask,
	GenVarGenBranchConstrMask                = 0x2000000000000    | GenericBranchingConstrMask,
            
	LimMemoryRankOneCutConstrMask            = 0x4000000000000    | InstMasterConstrMask | Base4NonLinearConstraintMask,
    //InstSeparationSpConstrMask               = 0x8000000000000    | SeparationSpConstrMask | InstanciatedConstrMask,
    //SecondStageCostSepSpConstrMask           = 0x10000000000000   | SeparationSpConstrMask | ConstraintMask,
    //SecondStageCostMastVarMask               = 0x20000000000000   | MasterVarMask | VariableMask,
    ExtendedArcCutConstrMask                 = 0x40000000000000   | InstMasterConstrMask | Base4NonLinearConstraintMask,
    CustomNonLinearCutConstrMask             = 0x80000000000000   | InstMasterConstrMask | Base4NonLinearConstraintMask,
    //BendersCutMask                           = 0x100000000000000  | InstMasterConstrMask,
	ResConsKnapsackCutConstrMask             = 0x200000000000000  | InstMasterConstrMask | Base4NonLinearConstraintMask,
	SoftConflictsCutConstrMask               = 0x400000000000000  | InstMasterConstrMask | Base4NonLinearConstraintMask,
	PackSetResConsInstMastBranchConstrMask   = 0x800000000000000  | InstMasterBranchingConstrMask
											                      | Base4NonLinearConstraintMask,
	CliqueCutConstrMask                      = 0x1000000000000000 | InstMasterConstrMask | Base4NonLinearConstraintMask,
	PackSetRyanFostInstMastBranchConstrMask  = 0x2000000000000000 | InstMasterBranchingConstrMask
												                  | Base4NonLinearConstraintMask,
    LimMemoryKPathCutConstrMask              = 0x4000000000000000 | InstMasterConstrMask | Base4NonLinearConstraintMask
	};
};


struct PcId
{
#ifdef _WINDOWS
  enum PcIdentifier : int64_t
#else
  enum PcIdentifier
#endif
    {
      MasterMask                    = 0x01,
	  SubproblemMask                            = 0x02,
	  ColGenSpConfMask                            = 0x04  | SubproblemMask
	};
};



template<typename Identifier>
inline bool compareIdentifier(const Identifier& type1, const Identifier& type2)
{
#ifdef _WINDOWS
  return ((type1 & type2) == (int64_t) type2);
#else
  return ((type1 & type2) == type2);
#endif
}


#endif /* BCVCIDENTIFIERC_HPP_ */
