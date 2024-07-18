/**
 *
 * This file bcMultiIndexC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMULTIINDEXC_H_
#define BCMULTIINDEXC_H_

#include <string>
#include "bcBoundLevC.hpp"

struct BcVarCoef;
class MultiIndex;
struct BcVarIndex;

struct BcVarConstrType
{
	enum BcVcType
	{
		local2Formulation = 0, globalSumOfLocals = 1, globalCloneOfLocals = 2
	};
};


class MultiIndexNames
{
  friend class MultiIndex;
  static const int maxPosition=8;
  char letterArray[maxPosition];
 public:
  MultiIndexNames(const char & firstName = '_',
                  const char & secondName = '_',
                  const char & thirdName = '_',
                  const char & fourthName = '_',
                  const char & fifthName =  '_',
                  const char & sixthName = '_',
                  const char & seventhName = '_',
                  const char & eighthName = '_');
  void append(const MultiIndexNames &idx);// added by Boris
  char first() const;
  char second() const;
  char third() const;
  char fourth() const;
  char fifth() const;
  char sixth() const;
  char seventh() const;
  char eighth() const;
};

class MultiIndex
{
	static const int maxPosition=8; 
	friend struct BcVarIndex;
	int indexArray[maxPosition];
public:
	int endPosition;
	MultiIndex(const MultiIndex & id);
	MultiIndex(const MultiIndex & id, int numIndicesToDelete);
	MultiIndex();
	MultiIndex(int first);
	MultiIndex(int first,int second);
	MultiIndex(int first,int second,int third);
	MultiIndex(int first,int second,int third,int fourth);
	MultiIndex(int first,int second,int third,int fourth,int fifth);
	MultiIndex(int first,int second,int third,int fourth,int fifth,int sixth);
	MultiIndex(int first,int second,int third,int fourth,int fifth,int sixth, int seventh);
	MultiIndex(int first,int second,int third,int fourth,int fifth,int sixth, int seventh, int eighth);
	void resetEnd();
	int index(int position) const;
	int first() const;
	int second() const;
	int third() const;
	int fourth() const;
	int fifth() const;
	int sixth() const;
	int seventh() const;
	int eighth() const;
	MultiIndex & operator+=(int index);
	bool operator<(const MultiIndex & that) const;
	bool operator==(const MultiIndex & that) const;
	const std::string & appendRef2name(std::string & name,
									   const MultiIndexNames & indexNames = MultiIndexNames()) const;
	const std::string & appendRefWithBrackets2name(std::string & name) const;
	friend std::ostream& operator<<(std::ostream& os, const MultiIndex & that);

    std::string indicesAsString() const; 
    
	// added by Boris
	/// \brief compute next mlti-index within the given range
	/// @return false iff this is equal to maxIndex 
	bool moveToNext(const MultiIndex &minIndex,const MultiIndex &maxIndex);
	/// \brief appends a set of indices to a multiindex
	void append(const MultiIndex& that);
	/// \returns the concatenation of this and that
	MultiIndex operator+(const MultiIndex &that) const;
};
std::ostream& operator<<(std::ostream& os, const MultiIndex & that);




#endif /* BCMULTIINDEXC_H_ */
